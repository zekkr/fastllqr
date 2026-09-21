! Weight production is separate from solver decisions; no LLQR/TVCQR flag.
subroutine ssqr_kernel_path(a,y,n,q,grid,ne,coordinate,kernel,h,tau,tol,maxit,threshold,min_keep,cache_flags,provider_flags, &
                            beta,hseq,diagnostics,ierr,failed_eval)
  use weighted_qr_core, only: fit_screen_path
  use iso_fortran_env, only: real64,int64
  implicit none
  integer,intent(in)::n,q,ne,kernel,maxit,min_keep,cache_flags,provider_flags
  real(real64),intent(in)::a(n,q),y(n),grid(ne),coordinate(n),h,tau,tol,threshold(ne)
  real(real64),intent(out)::beta(ne,q)
  integer,intent(out)::hseq(ne,q),diagnostics(ne,18),ierr,failed_eval
  real(real64),allocatable::sorted_coordinate(:),sort_work(:)
  integer,allocatable::sorted_ids(:),id_work(:),previous_active(:)
  integer(int64),allocatable::active_bits(:)
  integer::previous_na
  logical::coordinate_ordered,use_window
  if(h<=0 .or. (kernel/=1 .and. kernel/=2) .or. provider_flags<0 .or. provider_flags>3) then
    ierr=5; failed_eval=0; beta=0; hseq=0; diagnostics=0; return
  end if
  coordinate_ordered=.true.
  do previous_na=2,n
    if(coordinate(previous_na)<coordinate(previous_na-1)) coordinate_ordered=.false.
  end do
  previous_na=0
  use_window=kernel==2 .and. btest(provider_flags,0) .and. &
             (coordinate_ordered .or. btest(provider_flags,1))
  if(use_window) then
    allocate(sorted_coordinate(n),sort_work(n),sorted_ids(n),id_work(n),previous_active(n),active_bits((n+63)/64))
    if(coordinate_ordered) then
      sorted_coordinate=coordinate
      do previous_na=1,n
        sorted_ids(previous_na)=previous_na
      end do
      previous_na=0
    else
      call sort_coordinates()
    end if
    active_bits=0_int64
  end if
  call fit_screen_path(a,y,n,q,ne,tau,tol,maxit,threshold,min_keep,cache_flags,weights,beta,hseq,diagnostics,ierr,failed_eval)
contains
  subroutine weights(ev,nw,w,ids,na)
    integer,intent(in)::ev,nw
    real(real64),intent(inout)::w(nw)
    integer,intent(out)::ids(nw),na
    integer::i,ii,lo,hi,word,bitpos
    integer(int64)::bitword
    real(real64)::u,guard,left,right
    if(use_window) then
      do ii=1,previous_na
        w(previous_active(ii))=0
      end do
      active_bits=0_int64
      guard=8*epsilon(1.0_real64)*(abs(grid(ev))+abs(h)+1)
      left=grid(ev)-h-guard;right=grid(ev)+h+guard
      lo=lower_rank(left);hi=upper_rank(right)
      na=0
      if(coordinate_ordered .and. lo<=hi) then
        do i=lo,hi
          u=(grid(ev)-coordinate(i))/h;w(i)=0
          if(abs(u)<=1) w(i)=0.75_real64*(1-u*u)
          if(w(i)>0) then
            na=na+1;ids(na)=i
          end if
        end do
      else if(lo<=hi) then
        do ii=lo,hi
          i=sorted_ids(ii);u=(grid(ev)-coordinate(i))/h
          w(i)=0
          if(abs(u)<=1) w(i)=0.75_real64*(1-u*u)
          if(w(i)>0) active_bits((i-1)/64+1)=ibset(active_bits((i-1)/64+1),mod(i-1,64))
        end do
      end if
      if(.not.coordinate_ordered) then
      na=0
      do word=1,size(active_bits)
        bitword=active_bits(word)
        do while(bitword/=0_int64)
          bitpos=trailz(bitword);i=(word-1)*64+bitpos+1
          if(i<=nw) then
            na=na+1;ids(na)=i
          end if
          bitword=ibclr(bitword,bitpos)
        end do
      end do
      end if
      previous_na=na
      if(na>0) previous_active(1:na)=ids(1:na)
      return
    end if
    na=0
    do i=1,nw
      u=(grid(ev)-coordinate(i))/h
      if(kernel==1) then
        w(i)=exp(-0.5_real64*u*u)/sqrt(2*acos(-1.0_real64))
      else
        w(i)=0
        if(abs(u)<=1) w(i)=0.75_real64*(1-u*u)
      end if
      if(w(i)>0) then
        na=na+1; ids(na)=i
      end if
    end do
  end subroutine
  subroutine sort_coordinates()
    integer::i,j,k,left_start,mid,right_end,width
    sorted_coordinate=coordinate
    do i=1,n
      sorted_ids(i)=i
    end do
    width=1
    do while(width<n)
      left_start=1
      do while(left_start<=n)
        mid=min(left_start+width-1,n);right_end=min(left_start+2*width-1,n)
        i=left_start;j=mid+1;k=left_start
        do while(i<=mid .and. j<=right_end)
          if(sorted_coordinate(i)<=sorted_coordinate(j)) then
            sort_work(k)=sorted_coordinate(i);id_work(k)=sorted_ids(i);i=i+1
          else
            sort_work(k)=sorted_coordinate(j);id_work(k)=sorted_ids(j);j=j+1
          end if
          k=k+1
        end do
        do while(i<=mid)
          sort_work(k)=sorted_coordinate(i);id_work(k)=sorted_ids(i);i=i+1;k=k+1
        end do
        do while(j<=right_end)
          sort_work(k)=sorted_coordinate(j);id_work(k)=sorted_ids(j);j=j+1;k=k+1
        end do
        left_start=right_end+1
      end do
      sorted_coordinate=sort_work;sorted_ids=id_work;width=width*2
    end do
  end subroutine
  integer function lower_rank(value)
    real(real64),intent(in)::value
    integer::lo0,hi0,mid0
    lo0=1;hi0=n+1
    do while(lo0<hi0)
      mid0=lo0+(hi0-lo0)/2
      if(mid0<=n .and. sorted_coordinate(mid0)<value) then
        lo0=mid0+1
      else
        hi0=mid0
      end if
    end do
    lower_rank=lo0
  end function
  integer function upper_rank(value)
    real(real64),intent(in)::value
    integer::lo0,hi0,mid0
    lo0=0;hi0=n
    do while(lo0<hi0)
      mid0=lo0+(hi0-lo0+1)/2
      if(sorted_coordinate(mid0)<=value) then
        lo0=mid0
      else
        hi0=mid0-1
      end if
    end do
    upper_rank=lo0
  end function
end subroutine

! Arbitrary supplied weights expose the same core for model-independent tests.
subroutine ssqr_weighted_path(a,y,n,q,ne,weight_matrix,tau,tol,maxit,threshold,min_keep,cache_flags, &
                              beta,hseq,diagnostics,ierr,failed_eval)
  use weighted_qr_core, only: fit_screen_path
  use iso_fortran_env, only: real64
  implicit none
  integer,intent(in)::n,q,ne,maxit,min_keep,cache_flags
  real(real64),intent(in)::a(n,q),y(n),weight_matrix(n,ne),tau,tol,threshold(ne)
  real(real64),intent(out)::beta(ne,q)
  integer,intent(out)::hseq(ne,q),diagnostics(ne,18),ierr,failed_eval
  call fit_screen_path(a,y,n,q,ne,tau,tol,maxit,threshold,min_keep,cache_flags,weights,beta,hseq,diagnostics,ierr,failed_eval)
contains
  subroutine weights(ev,nw,w,ids,na)
    integer,intent(in)::ev,nw
    real(real64),intent(inout)::w(nw)
    integer,intent(out)::ids(nw),na
    integer::i
    w=weight_matrix(:,ev); na=0
    do i=1,nw
      if(w(i)>0) then
        na=na+1; ids(na)=i
      end if
    end do
  end subroutine
end subroutine
