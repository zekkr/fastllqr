! Common seq-screen weighted QR; derived from LLQR v64 / TVCQR v41.
! Only fixed A, y, a weight provider, thresholds and numerical controls enter.
module weighted_qr_core
  use iso_fortran_env, only: real64,int64
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: fit_screen_path
  integer, parameter :: dp=real64
  abstract interface
    subroutine weight_provider(ev,n,w,ids,na)
      import dp
      integer,intent(in)::ev,n
      real(dp),intent(inout)::w(n)
      integer,intent(out)::ids(n),na
    end subroutine
  end interface
contains
  subroutine fit_screen_path(a,y,n,q,ne,tau,tol,maxit,threshold,min_keep,cache_flags,get_weights, &
                             beta_out,h_out,diagnostics,ierr,failed_eval)
    integer,intent(in)::n,q,ne,maxit,min_keep,cache_flags
    real(dp),intent(in)::a(n,q),y(n),tau,tol,threshold(ne)
    procedure(weight_provider)::get_weights
    real(dp),intent(out)::beta_out(ne,q)
    integer,intent(out)::h_out(ne,q),diagnostics(ne,18),ierr,failed_eval
    real(dp),allocatable::gx(:,:),bv(:),ws(:),rx(:,:),ry(:),w(:),rprev(:),rcand(:)
    real(dp),allocatable::inverse_work(:)
    real(dp),allocatable::raw_row_cache(:,:),raw_rhs_cache(:),raw_sign_cache(:)
    real(dp)::inv(q,q),beta(q),previous_beta(q),xh(q,q),aggx(2,q),aggy(2),aggw(2)
    real(dp)::inverse_cache(q,q)
    real(dp)::gamma,pivot_tol,res_tol,work_query(1),row_sign,rr,largest
    integer::hp(q),hc(q),hsub(q),r1(q),r2(q),inv_piv(q),lwork,info
    integer::inverse_owner(q),row_owner(q)
    integer::active_ids(n),kept(n),map(n),tag(n),bad(n),ib(n+3)
    integer::ev,na,ns,ms,i,j,k,ii,g,nb,attempt,iters,cold_status,solve_status,nh
    integer(int64)::residual_epoch,residual_stamp(n),previous_bits(q),candidate_bits(q)
    integer(int64)::row_epoch,row_stamp(n)
    logical::fvr(n+3),in_h(n),reclassify,force_full,ok,warm_ok,done,hasagg(2)
    logical::reuse_residuals,use_certificate,audit_certificate,same_image,certificate_hit
    logical::reuse_raw_rows,full_recorded,use_active_lists

    ierr=0; failed_eval=0; beta_out=0; h_out=0; diagnostics=0
    diagnostics(:,13:18)=-1
    if(n<q .or. q<1 .or. ne<1 .or. min_keep<1 .or. maxit<1) then
      ierr=5; return
    end if
    if(tau<=0 .or. tau>=1 .or. tol<=0 .or. any(threshold<=0)) then
      ierr=5; return
    end if
    if(.not.all(ieee_is_finite(a)) .or. .not.all(ieee_is_finite(y)) .or. &
       .not.all(ieee_is_finite(threshold))) then
      ierr=5; return
    end if
    allocate(gx(n+3,q),bv(n+3),ws(n+2),rx(n+2,q),ry(n+2),w(n),rprev(n),rcand(n))
    ! Fit-wide inverse workspace (TVCQR v41); q=2 keeps v64 arithmetic.
    xh=0; inv_piv=1
    call dgetri(q,xh,q,inv_piv,work_query,-1,info)
    if(info/=0) then
      ierr=5; return
    end if
    lwork=max(q,int(work_query(1))); allocate(inverse_work(lwork))
    pivot_tol=max(10*tol,1.0e-12_dp); res_tol=1.0e-6_dp
    previous_beta=0; hp=0;in_h=.false.
    reuse_residuals=btest(cache_flags,0);use_certificate=btest(cache_flags,1)
    audit_certificate=btest(cache_flags,2)
    reuse_raw_rows=btest(cache_flags,3)
    use_active_lists=btest(cache_flags,4)
    if(cache_flags<0 .or. cache_flags>31 .or. (use_certificate .and. .not.reuse_residuals)) then
      ierr=5; return
    end if
    residual_stamp=0_int64;residual_epoch=1_int64
    w=0
    inverse_owner=0;row_owner=0;row_epoch=0_int64;row_stamp=0_int64
    if(reuse_raw_rows) allocate(raw_row_cache(n,q),raw_rhs_cache(n),raw_sign_cache(n))

    do ev=1,ne
      call get_weights(ev,n,w,active_ids,na)
      diagnostics(ev,1)=na
      if(na<q .or. any(w<0) .or. .not.all(ieee_is_finite(w))) then
        call fail(3); return
      end if
      if(ev==1) then
        ! Shared active-first objective. The more expensive guarded full-n
        ! interpolating basis is a recovery option, not a per-model fast path.
        ns=na; kept(1:ns)=active_ids(1:na)
        call pack_individuals()
        ms=ns
        call cold_tableau(previous_beta,0,solve_status)
        if(solve_status/=0) then
          cold_status=solve_status;diagnostics(ev,7)=4
          call seeded_first(solve_status)
          if(solve_status==0) then
            diagnostics(ev,6)=4;diagnostics(ev,12)=1
          else if(na<n) then
          diagnostics(ev,12)=1
          ns=n
          do i=1,n
            kept(i)=i
          end do
          call pack_individuals()
          ms=ns
          call cold_tableau(previous_beta,0,solve_status)
          else
            solve_status=cold_status
          end if
        end if
        if(solve_status/=0) then
          call fail(solve_status); return
        end if
        diagnostics(ev,2)=ms; diagnostics(ev,3)=ms
        beta_out(ev,:)=beta; h_out(ev,:)=hc
        previous_beta=beta; hp=hc
        previous_bits=transfer(beta,previous_bits)
        cycle
      end if

      ! Materialize only rows consumed at this transition. Never use stale rows.
      do ii=1,na
        i=active_ids(ii)
        if(.not.reuse_residuals .or. residual_stamp(i)/=residual_epoch) then
          rprev(i)=raw_residual(i,previous_beta)
          residual_stamp(i)=residual_epoch
          diagnostics(ev,11)=diagnostics(ev,11)+1
        end if
        if(.not.ieee_is_finite(rprev(i))) then
          call fail(4); return
        end if
      end do
      call invert_h(hp,inv,ok)
      if(.not.use_active_lists) in_h=.false.
      if(ok) in_h(hp)=.true.
      warm_ok=ok
      if(.not.warm_ok) then
        diagnostics(ev,6)=2; diagnostics(ev,7)=1
      end if
      gamma=threshold(ev); reclassify=.true.; force_full=.false.; done=.false.
      full_recorded=.false.
      attempt=0
      do while(.not.done)
        attempt=attempt+1
        if(attempt>max(10,2*q)) force_full=.true.
        if(force_full .and. .not.full_recorded) then
          ! Count abandonment of screening even when the warm full-active solve
          ! succeeds. This is recovery, not a reduced solve and not seq fallback.
          diagnostics(ev,6)=3;diagnostics(ev,8)=diagnostics(ev,8)+1
          full_recorded=.true.
        end if
        if(reclassify .or. force_full) then
          if(use_active_lists) then
            tag(active_ids(1:na))=0
          else
            tag=0
          end if
          if(.not.force_full) then
            nh=0; largest=0
            do ii=1,na
              i=active_ids(ii)
              if(abs(rprev(i))<=gamma) nh=nh+1
              largest=max(largest,abs(rprev(i)))
            end do
            if(attempt==1) then
              diagnostics(ev,13)=nh
              diagnostics(ev,14)=0
            end if
            do while(nh<min(min_keep,na) .and. gamma<largest)
              gamma=gamma*1.5_dp; nh=0
              if(attempt==1) diagnostics(ev,14)=diagnostics(ev,14)+1
              do ii=1,na
                if(abs(rprev(active_ids(ii)))<=gamma) nh=nh+1
              end do
            end do
            if(attempt==1) diagnostics(ev,15)=nh
            do ii=1,na
              i=active_ids(ii)
              if(rprev(i)<-gamma) tag(i)=-1
              if(rprev(i)>gamma) tag(i)=1
            end do
          end if
        end if
        if(warm_ok) tag(hp)=0
        ! Ordered individual packing and exact current-weight aggregation.
        ns=0; aggx=0; aggy=0; aggw=0
        do ii=1,merge(na,n,use_active_lists)
          i=ii
          if(use_active_lists) i=active_ids(ii)
          if(w(i)<=0 .and. .not.in_h(i)) cycle
          if(tag(i)==0) then
            ns=ns+1; kept(ns)=i
          else
            g=1
            if(tag(i)>0) g=2
            aggw(g)=aggw(g)+w(i)
            do j=1,q
              aggx(g,j)=aggx(g,j)+a(i,j)*w(i)
            end do
            aggy(g)=aggy(g)+y(i)*w(i)
          end if
        end do
        if(use_active_lists .and. warm_ok) then
          ! Providers return ascending original IDs. Insert zero-weight previous
          ! H padding in that order without traversing the other inactive rows.
          do j=1,q
            i=hp(j)
            if(w(i)>0) cycle
            k=ns
            do while(k>0)
              if(kept(k)<i) exit
              kept(k+1)=kept(k);k=k-1
            end do
            kept(k+1)=i;ns=ns+1
          end do
        end if
        call pack_individuals()
        ms=ns; hasagg=aggw>0
        do g=1,2
          if(.not.hasagg(g)) cycle
          ms=ms+1; rx(ms,:)=aggx(g,:); ry(ms)=aggy(g); ws(ms)=1
        end do
        if(attempt==1) then
          diagnostics(ev,2)=ms
          diagnostics(ev,16)=ns-diagnostics(ev,15)
          diagnostics(ev,17)=count(hasagg)
          diagnostics(ev,18)=ns
        end if
        diagnostics(ev,3)=ms
        if(ns<q) then
          solve_status=3
        else
          if(warm_ok) then
            call fresh_tableau(ok)
          else
            ok=.false.
          end if
          if(ok) then
            call simplex(gx,bv,ib,fvr,r1,r2,ws,n+3,ms,q,tol,pivot_tol,maxit,iters,solve_status)
            diagnostics(ev,4)=diagnostics(ev,4)+iters
            if(solve_status==0) call refit_candidate(solve_status)
          else
            if(.not.force_full) diagnostics(ev,6)=2
            diagnostics(ev,9)=diagnostics(ev,9)+1
            call cold_tableau(previous_beta,ns,solve_status)
          end if
        end if
        if(solve_status/=0) then
          if(force_full) then
            ! Explicit ppro full-active cold recovery, not a call to seq.
            ns=na; kept(1:ns)=active_ids(1:na); call pack_individuals(); ms=ns
            diagnostics(ev,6)=3
            call cold_tableau(previous_beta,0,solve_status)
            if(solve_status/=0) then
              call fail(solve_status); return
            end if
            tag=0; diagnostics(ev,3)=ms
          else
            diagnostics(ev,5)=diagnostics(ev,5)+1
            gamma=gamma*2; reclassify=.true.
            if(attempt>=3) force_full=.true.
            cycle
          end if
        end if
        candidate_bits=transfer(beta,candidate_bits)
        same_image=all(candidate_bits==previous_bits)
        if(same_image) then
          do i=1,q
            if(.not.any(hc(i)==hp)) same_image=.false.
          end do
        end if
        ! Common exact-image rule: no approximate beta comparison. All current
        ! rows were materialized above from this accepted owner. H membership
        ! is checked conservatively, although these buffers keep raw H residuals.
        certificate_hit=use_certificate .and. same_image .and. ieee_is_finite(gamma) .and. gamma>=0
        nb=0
        if(certificate_hit) then
          diagnostics(ev,10)=diagnostics(ev,10)+1
          if(audit_certificate) then
            do ii=1,na
              i=active_ids(ii);rr=raw_residual(i,beta)
              if(transfer(rr,0_int64)/=transfer(rprev(i),0_int64)) error stop 'exact residual image differs'
              if((tag(i)==1 .and. rr<=0) .or. (tag(i)==-1 .and. rr>=0)) &
                error stop 'exact residual certificate contradicts literal signs'
            end do
          end if
        else
        ! Literal strict omitted-row verification. Exact zero is a bad sign.
        do ii=1,na
          i=active_ids(ii)
          rcand(i)=raw_residual(i,beta)
          if((tag(i)==1 .and. rcand(i)<=0) .or. (tag(i)==-1 .and. rcand(i)>=0)) then
            nb=nb+1; bad(nb)=i
          end if
        end do
        diagnostics(ev,11)=diagnostics(ev,11)+na
        if(.not.all(ieee_is_finite(rcand(active_ids(1:na))))) then
          call fail(4); return
        end if
        end if
        if(nb>0) then
          diagnostics(ev,5)=diagnostics(ev,5)+1
          if(real(nb,dp)>0.1_dp*real(ms,dp)) then
            gamma=gamma*2; reclassify=.true.
          else
            tag(bad(1:nb))=0; reclassify=.false.
          end if
        else
          done=.true.
        end if
      end do
      beta_out(ev,:)=beta; h_out(ev,:)=hc
      if(reuse_residuals .and. .not.same_image) then
        residual_epoch=residual_epoch+1_int64
        do ii=1,na
          i=active_ids(ii)
          rprev(i)=rcand(i);residual_stamp(i)=residual_epoch
        end do
      end if
      if(warm_ok) in_h(hp)=.false.
      previous_beta=beta; hp=hc
      previous_bits=candidate_bits
    end do
  contains
    subroutine fail(code)
      integer,intent(in)::code
      ierr=code; failed_eval=ev
    end subroutine
    real(dp) function raw_residual(id,b)
      integer,intent(in)::id
      real(dp),intent(in)::b(q)
      integer::col
      raw_residual=y(id)
      if(q==2) then
        raw_residual=raw_residual-a(id,1)*b(1)
        raw_residual=raw_residual-a(id,2)*b(2)
        return
      end if
      do col=1,q
        raw_residual=raw_residual-a(id,col)*b(col)
      end do
    end function
    subroutine invert_h(hh,ainv,success)
      integer,intent(in)::hh(q)
      real(dp),intent(out)::ainv(q,q)
      logical,intent(out)::success
      integer::ci,cj,inf,ip(q)
      real(dp)::det
      success=.false.
      if(reuse_raw_rows .and. all(hh==inverse_owner)) then
        ainv=inverse_cache;success=.true.;return
      end if
      do ci=1,q
        if(hh(ci)<1 .or. hh(ci)>n) return
        do cj=ci+1,q
          if(hh(ci)==hh(cj)) return
        end do
        ainv(ci,:)=a(hh(ci),:)
      end do
      if(q==2) then
        det=ainv(1,1)*ainv(2,2)-ainv(1,2)*ainv(2,1)
        if(abs(det)<=1.0e-10_dp) return
        xh=ainv
        ainv(1,1)=xh(2,2)/det; ainv(1,2)=-xh(1,2)/det
        ainv(2,1)=-xh(2,1)/det; ainv(2,2)=xh(1,1)/det
      else
        call dgetrf(q,q,ainv,q,ip,inf)
        if(inf/=0) return
        call dgetri(q,ainv,q,ip,inverse_work,lwork,inf)
        if(inf/=0) return
      end if
      success=all(ieee_is_finite(ainv))
      if(success .and. reuse_raw_rows) then
        inverse_owner=hh;inverse_cache=ainv
      end if
    end subroutine
    subroutine pack_individuals()
      integer::ci
      if(use_active_lists) then
        do ci=1,q
          if(hp(ci)>=1 .and. hp(ci)<=n) map(hp(ci))=0
        end do
      else
        map=0
      end if
      do ci=1,ns
        map(kept(ci))=ci
        rx(ci,:)=a(kept(ci),:); ry(ci)=y(kept(ci)); ws(ci)=w(kept(ci))
      end do
    end subroutine
    subroutine seeded_first(status)
      integer,intent(out)::status
      integer::seed_h(q),ci,cj,ck,pass,best_id,row
      real(dp)::orth(q,q),vrow(q),vbest(q),normrow,normbest,best_weight,max_weight,projection
      real(dp)::first_inverse(q,q),first_beta(q),tr(q),rhs,sgn,cost
      logical::valid
      status=3;orth=0;seed_h=0;max_weight=maxval(w)
      ! Weighted, rank-revealing row selection. Two projection passes guard
      ! nearly dependent rows; the same rule is used at every q and model.
      do cj=1,q
        best_id=0;best_weight=-1;normbest=0
        do ci=1,n
          if(w(ci)<=best_weight .or. any(seed_h==ci)) cycle
          vrow=a(ci,:)
          do pass=1,2
            do ck=1,cj-1
              projection=sum(vrow*orth(ck,:))
              vrow=vrow-projection*orth(ck,:)
            end do
          end do
          normrow=sqrt(sum(vrow*vrow))
          if(normrow<=1.0e-10_dp*max(1.0_dp,sqrt(sum(a(ci,:)**2)))) cycle
          best_id=ci;best_weight=w(ci);vbest=vrow;normbest=normrow
        end do
        if(best_id==0) return
        if(best_weight<=max(sqrt(epsilon(1.0_dp))*max_weight,128*tol)) return
        seed_h(cj)=best_id;orth(cj,:)=vbest/normbest
      end do
      call invert_h(seed_h,first_inverse,valid)
      if(.not.valid) return
      first_beta=0
      do ci=1,q
        do cj=1,q
          first_beta(ci)=first_beta(ci)+first_inverse(ci,cj)*y(seed_h(cj))
        end do
      end do
      if(.not.all(ieee_is_finite(first_beta))) return
      ns=n;ms=n
      do ci=1,n
        kept(ci)=ci
      end do
      call pack_individuals()
      gx(1:q,:)=first_inverse;bv(1:q)=first_beta;fvr(1:ms+1)=.false.
      fvr(1:q)=.true.;fvr(ms+1)=.true.;ib(ms+1)=0;bv(ms+1)=0
      do ci=1,q
        ib(ci)=ci;r1(ci)=q+seed_h(ci);r2(ci)=r1(ci)+ms
        gx(ms+1,ci)=tau*w(seed_h(ci))
      end do
      row=q
      do ci=1,n
        if(any(seed_h==ci)) cycle
        row=row+1;rhs=raw_residual(ci,first_beta);sgn=1
        if(rhs<0) sgn=-1
        tr=0
        do cj=1,q
          do ck=1,q
            tr(cj)=tr(cj)+a(ci,ck)*first_inverse(ck,cj)
          end do
        end do
        gx(row,:)=-sgn*tr;bv(row)=sgn*rhs
        ib(row)=q+ci;cost=tau*w(ci)
        if(sgn<0) then
          ib(row)=q+ms+ci;cost=(1-tau)*w(ci)
        end if
        gx(ms+1,:)=gx(ms+1,:)-cost*gx(row,:)
      end do
      call simplex(gx,bv,ib,fvr,r1,r2,ws,n+3,ms,q,tol,pivot_tol,maxit,iters,status)
      diagnostics(ev,4)=diagnostics(ev,4)+iters
      if(status==0) call refit_candidate(status)
    end subroutine
    subroutine refit_candidate(status)
      integer,intent(out)::status
      integer::ci,cj,hi
      logical::valid
      real(dp)::ci_inv(q,q)
      status=3
      do ci=1,q
        hi=r1(ci)-q
        if(hi<1 .or. hi>ns) return
        hc(ci)=kept(hi)
      end do
      call invert_h(hc,ci_inv,valid)
      if(.not.valid) return
      beta=0
      do ci=1,q
        do cj=1,q
          beta(ci)=beta(ci)+ci_inv(ci,cj)*y(hc(cj))
        end do
      end do
      if(.not.all(ieee_is_finite(beta))) return
      do ci=1,q
        if(abs(raw_residual(hc(ci),beta))>res_tol) return
      end do
      status=0
    end subroutine
    subroutine cold_tableau(offset,n_individual,status)
      real(dp),intent(in)::offset(q)
      integer,intent(in)::n_individual
      integer,intent(out)::status
      integer::ci,cj
      real(dp)::rhs,cost
      gx(1:ms+1,:)=0; bv(1:ms+1)=0; fvr(1:ms+1)=.false.; ib(1:ms+1)=0
      do ci=1,ms
        rhs=ry(ci)
        do cj=1,q
          rhs=rhs-rx(ci,cj)*offset(cj)
        end do
        row_sign=1
        if(rhs<0) row_sign=-1
        if(n_individual>0 .and. ci>n_individual) then
          fvr(ci)=.true.
          row_sign=1
          if(hasagg(1) .and. ci==n_individual+1) row_sign=-1
        end if
        gx(ci,:)=row_sign*rx(ci,:); bv(ci)=row_sign*rhs
        ib(ci)=q+ci; cost=tau*ws(ci)
        if(row_sign<0) then
          ib(ci)=q+ms+ci; cost=(1-tau)*ws(ci)
        end if
        gx(ms+1,:)=gx(ms+1,:)-cost*gx(ci,:)
      end do
      fvr(ms+1)=.true.
      do ci=1,q
        r1(ci)=ci
      end do
      r2=0
      call simplex(gx,bv,ib,fvr,r1,r2,ws,n+3,ms,q,tol,pivot_tol,maxit,iters,status)
      diagnostics(ev,4)=diagnostics(ev,4)+iters
      if(status==0) call refit_candidate(status)
    end subroutine
    subroutine fresh_tableau(success)
      logical,intent(out)::success
      integer::ci,cj,ck,row,pass,id,source_row,first_source,last_source
      real(dp)::tr(q),rhs,lambda
      logical::row_hit
      success=.false.
      if(reuse_raw_rows .and. any(hp/=row_owner)) then
        row_owner=hp;row_epoch=row_epoch+1_int64
      end if
      do ci=1,q
        hsub(ci)=map(hp(ci))
        if(hsub(ci)==0) return
      end do
      gx(1:q,:)=inv; bv(1:q)=0; fvr(1:ms+1)=.false.; ib(1:ms+1)=0
      do ci=1,q
        ib(ci)=ci; fvr(ci)=.true.
        r1(ci)=q+hsub(ci); r2(ci)=r1(ci)+ms
        do cj=1,q
          bv(ci)=bv(ci)+inv(ci,cj)*y(hp(cj))
        end do
        gx(ms+1,ci)=tau*ws(hsub(ci))
      end do
      bv(ms+1)=0; fvr(ms+1)=.true.; row=q
      ! Match the parents' positive, negative, low-aggregate, high-aggregate order.
      do pass=1,4
        if(pass<=2) then
          first_source=1;last_source=ns
        else if(pass==3 .and. hasagg(1)) then
          first_source=ns+1;last_source=ns+1
        else if(pass==4 .and. hasagg(2)) then
          first_source=ms;last_source=ms
        else
          cycle
        end if
        do source_row=first_source,last_source
          if(source_row<=ns) then
            id=kept(source_row)
            if(in_h(id)) cycle
            if(pass>2) cycle
            rr=rprev(id)
            if(pass==1 .and. rr<0) cycle
            if(pass==2 .and. rr>=0) cycle
            row_sign=1
            if(pass==2) row_sign=-1
          else
            if(pass<=2) cycle
            if(pass==3) then
              if(.not.hasagg(1) .or. source_row/=ns+1) cycle
              row_sign=-1
            else
              if(.not.hasagg(2) .or. source_row/=ms) cycle
              row_sign=1
            end if
          end if
          row=row+1
          if(row>ms) return
          ! Only raw, pre-pivot individual rows with identical ordered H and
          ! sign are reusable. Current weighted aggregates are always rebuilt.
          row_hit=.false.
          if(reuse_raw_rows .and. source_row<=ns) then
            if(row_stamp(id)==row_epoch) row_hit=(raw_sign_cache(id)==row_sign)
          end if
          if(row_hit) then
            tr=raw_row_cache(id,:);rhs=raw_rhs_cache(id)
          else
            tr=0
            do cj=1,q
              do ck=1,q
                tr(cj)=tr(cj)+row_sign*rx(source_row,ck)*inv(ck,cj)
              end do
            end do
            rhs=0
            do cj=1,q
              rhs=rhs-tr(cj)*y(hp(cj))
            end do
            rhs=rhs+row_sign*ry(source_row)
            if(reuse_raw_rows .and. source_row<=ns) then
              raw_row_cache(id,:)=tr;raw_rhs_cache(id)=rhs
              raw_sign_cache(id)=row_sign;row_stamp(id)=row_epoch
            end if
          end if
          gx(row,:)=-tr; bv(row)=rhs
          ib(row)=q+source_row; lambda=tau*ws(source_row)
          if(row_sign<0) then
            ib(row)=q+ms+source_row; lambda=(1-tau)*ws(source_row)
          end if
          if(source_row>ns) fvr(row)=.true.
          do cj=1,q
            gx(ms+1,cj)=gx(ms+1,cj)+lambda*tr(cj)
          end do
        end do
      end do
      success=(row==ms)
    end subroutine
  end subroutine fit_screen_path

  ! One residual-pair simplex implementation for both cold and H-block starts.
  ! LLQR v64 algebra, generalized q; TVCQR v41 pivot safety policy shared by all q.
  subroutine simplex(gx,bv,ib,fvr,r1,r2,w,ld,ms,q,tol,ptol,maxit,iters,status)
    integer,intent(in)::ld,ms,q,maxit
    real(dp),intent(inout)::gx(ld,q),bv(ld)
    integer,intent(inout)::ib(ld),r1(q),r2(q)
    logical,intent(inout)::fvr(ld)
    real(dp),intent(in)::w(ms),tol,ptol
    integer,intent(out)::iters,status
    real(dp)::rr(2,q),yy(ms+1),ee(ms+1),minimum,ratio,minratio,pivot
    integer::i,j,k,col,side,enter,idx
    status=2; iters=0
    do while(iters<maxit)
      do j=1,q
        rr(1,j)=gx(ms+1,j)
        if(r2(j)==0) then
          rr(1,j)=-abs(rr(1,j)); rr(2,j)=0
        else
          idx=r1(j)-q
          if(idx<1 .or. idx>ms) then
            status=5; return
          end if
          rr(2,j)=w(idx)-rr(1,j)
        end if
      end do
      minimum=minval(rr)
      if(.not.ieee_is_finite(minimum)) then
        status=4; return
      end if
      if(minimum>=-tol) then
        status=0; return
      end if
      col=1; side=1
      select_enter: do j=1,q
        do i=1,2
          if(rr(i,j)==minimum) then
            col=j; side=i; exit select_enter
          end if
        end do
      end do select_enter
      enter=r1(col)
      if(side==2) enter=r2(col)
      yy=gx(1:ms+1,col)
      if(r2(col)/=0 .and. side==2) yy=-yy
      minratio=huge(1.0_dp); k=0
      do i=1,ms
        if(fvr(i)) cycle
        if(r2(col)==0 .and. yy(ms+1)>=0) then
          if(yy(i)>=-ptol) cycle
          ratio=-bv(i)/yy(i)
        else
          if(yy(i)<=ptol) cycle
          ratio=bv(i)/yy(i)
        end if
        if(ratio<minratio-tol) then
          minratio=ratio; k=i
        end if
      end do
      if(k==0) return
      if(r2(col)==0) then
        fvr(k)=.true.
      else if(side==2) then
        yy(ms+1)=yy(ms+1)+w(r1(col)-q)
      end if
      do i=1,ms+1
        ee(i)=yy(i)/yy(k)
      end do
      ee(k)=1-1/yy(k)
      gx(1:ms+1,col)=0
      if(ib(k)<=q+ms) then
        gx(k,col)=1; r1(col)=ib(k); r2(col)=ib(k)+ms
      else
        gx(k,col)=-1; gx(ms+1,col)=w(ib(k)-ms-q)
        r1(col)=ib(k)-ms; r2(col)=ib(k)
      end if
      do j=1,q
        pivot=gx(k,j)
        do i=1,ms+1
          gx(i,j)=gx(i,j)-ee(i)*pivot
        end do
      end do
      pivot=bv(k)
      do i=1,ms+1
        bv(i)=bv(i)-ee(i)*pivot
      end do
      ib(k)=enter; iters=iters+1
    end do
  end subroutine simplex
end module weighted_qr_core
