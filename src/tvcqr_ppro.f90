! Time-varying coefficient quantile regression - Sequential preprocessing
! Faithful translation of tvcqr_seq_ppro.R with focus on correctness for eva_t <= 4
!
subroutine tvcqr_seq_ppro_fortran(x, y, m, nvar, tau, h, h_factor, tol, maxit, &
                                   bland_int, Mm_factor, &
                                   theta_ll_est, beta_full_est, H_seq, &
                                   full_active_recovery_count, ierr, failed_eval, min_subsample_size_in, &
                                   always_same_h_refit_int, threshold_lower_bound_int, &
                                   threshold_scale_mode_int)

    implicit none

    ! Input arguments
    integer, intent(in) :: m, nvar, maxit, bland_int
    integer, intent(in) :: min_subsample_size_in, always_same_h_refit_int
    integer, intent(in) :: threshold_lower_bound_int, threshold_scale_mode_int
    double precision, intent(in) :: x(m, nvar), y(m), tau, tol, h_factor
    double precision, intent(in) :: Mm_factor
    double precision, intent(inout) :: h

    ! Output arguments
    double precision, intent(out) :: theta_ll_est(m, nvar+1)
    double precision, intent(out) :: beta_full_est(m, 2*(nvar+1))
    integer, intent(out) :: H_seq(m, 2*(nvar+1))
    integer, intent(out) :: full_active_recovery_count
    integer, intent(out) :: ierr
    integer, intent(out) :: failed_eval

    ! Local variables
    double precision :: x_norms(m), mm_thresh, mmm_thresh, M_threshold, threshold_scale
    logical :: sl(m), sh(m), not_jl_or_jh(m), active(m)
    logical :: has_sl_agg, has_sh_agg
    integer :: idx_not_jl_or_jh(m)
    double precision :: temp_check
    integer :: min_subsample_size, n_potential_S
    integer :: n_active, target_min
    double precision :: residual_scale, pivot_tol
    double precision :: abs_r(m)  ! NEW - for median of absolute residuals

    integer :: ms, ms_org
    double precision :: ws(m+3)

    double precision :: glob_wx(2*(nvar+1)), glob_wy
    double precision :: ghib_wx(2*(nvar+1)), ghib_wy
    double precision :: wsl(m), wsh(m)

    !double precision :: gammaxs_temp(m+3, 2*(nvar+1))
    !double precision :: bs_temp(m+3)
    !double precision :: gammaxs(m+3, 2*(nvar+1))
    !double precision :: bs(m+3)
    double precision :: xhinv(2*(nvar+1), 2*(nvar+1))
    double precision :: Pxhbarxhinv(m, 2*(nvar+1))
    double precision :: lambda(m+3)

    double precision :: best_pivot_value
    integer :: best_k
    double precision :: ratio

    ! Storage for reuse in eva_t >= 3
    double precision :: gammaxs_pos(m, 2*(nvar+1))
    double precision :: gammaxs_neg(m, 2*(nvar+1))
    double precision :: bs_pos(m), bs_neg(m)
    integer :: idx_Hbar_pos(m), idx_Hbar_neg(m)
    integer :: n_Hbar_pos, n_Hbar_neg

    double precision :: time_index(m), w(m)
    double precision :: cc(2*(1+nvar) + 2*m)
    !double precision :: A(m, 2*(nvar+1))
    !double precision :: gammax(m+1, 2*(nvar+1))
    double precision :: b(m+1)
    integer :: IB(m+1), IBs(m+3)
    logical :: freevarrow(m+3)
    integer :: r1(2*(nvar+1)), r2(2*(nvar+1))
    double precision :: rr(2, 2*(nvar+1))

    double precision :: yy(m+3), ee(m+3), k_vals(m+3)
    double precision :: u(m), v(m), estimate(2*(nvar+1))
    double precision :: estimate_refit(2*(nvar+1)), r_refit(m)
    double precision :: r(m), r_prev(m)
    double precision :: beta_prev(2*(nvar+1))
    double precision :: pivot_row(2*(nvar+1))

    integer :: i, j, k, t, eva_t, iter
    integer :: t_rr, tsep
    double precision :: rrl, min_k, temp_sum
    logical :: bland, always_same_h_refit, threshold_lower_bound
    double precision :: b_k_original
    logical :: not_optimal, not_new_sl_sh, same_h_ok, refit_used
    integer :: bad_signs, n_sl, n_sh, n_sure_signs
    integer :: idpos(m), idneg(m), n_idpos, n_idneg
    integer :: H_indices(2*(nvar+1)), Hbar_indices(m+2)
    integer :: u_in_IBs(m), v_in_IBs(m)
    integer :: idx

    ! Variables for eva_t >= 3
    integer :: idx_Hbar_pos2(m), idx_Hbar_neg2(m)
    logical :: matched_rows_pos(m), matched_rows_neg(m)
    integer :: rows_pos(m), rows_neg(m)
    integer :: n_matched_pos, n_matched_neg
    integer :: id_gammaxs_Hbar(m)
    double precision :: temp_vec(2*(nvar+1))

    integer :: preprocessing_attempts
    integer :: inv_info
    logical :: valid_H
    double precision :: temp_vec1(2*(nvar+1))
    integer :: ii, kk
    logical :: force_full_sample, accept_subsample, no_pivot_flag
    logical :: use_independent_init, previous_H_valid, shifted_init_success
    logical :: full_recovery_success
    integer :: independent_trigger
    integer :: empty_pivot_count, max_empty_pivot_retries
    double precision :: res_tol


    integer :: n_hbar_rows        ! For the critical lambda/Pxhbarxhinv dimension
    integer :: idx_count, idx_loop
    double precision :: weight_sum, sum_w_sl, sum_w_sh

    double precision, allocatable :: A(:,:)           ! Size: m × 2(nvar+1)
    double precision, allocatable :: gammax(:,:)      ! Size: (m+1) × 2(nvar+1)
    double precision, allocatable :: gammaxs_temp(:,:)! Size: (m+3) × 2(nvar+1)
    double precision, allocatable :: bs_temp(:)       ! Size: m+3
    double precision, allocatable :: gammaxs(:,:)     ! Size: (m+3) × 2(nvar+1)
    double precision, allocatable :: bs(:)            ! Size: m+3
    logical :: unbounded_detected
    logical :: simplex_converged
    integer :: iter_attempt
    integer :: total_simplex_iterations
    integer :: total_preprocessing_loops
    integer :: max_iter_at_any_t





    ierr = 0
    failed_eval = 0
    total_simplex_iterations = 0
    total_preprocessing_loops = 0
    max_iter_at_any_t = 0




    ! Start of executable code

    ! Convert integer to logical for bland
    bland = (bland_int /= 0)
    always_same_h_refit = (always_same_h_refit_int /= 0)
    threshold_lower_bound = (threshold_lower_bound_int /= 0)
    pivot_tol = max(10.0d0 * tol, 1.0d-12)
    max_empty_pivot_retries = 3
    res_tol = 1.0d-6

    ! Allocate the big arrays
    allocate(A(m, 2*(nvar+1)))
    allocate(gammax(m+1, 2*(nvar+1)))
    allocate(gammaxs_temp(m+3, 2*(nvar+1)))
    allocate(bs_temp(m+3))
    allocate(gammaxs(m+3, 2*(nvar+1)))
    allocate(bs(m+3))

    ! Set default bandwidth if h = 0
    if (h <= 0.0d0) then
        h = dble(m)**(-0.2d0) * h_factor
    end if

    ! Calculate x norms
    do i = 1, m
        x_norms(i) = 0.0d0
        do j = 1, nvar
            x_norms(i) = x_norms(i) + x(i,j)**2
        end do
        x_norms(i) = sqrt(x_norms(i))
    end do

    ! Calculate mm threshold - matching R: mm <- log(m)^{4} * h^2 *max(x.norms)
    mm_thresh = log(dble(m))**4 * h**2 * maxval(x_norms)
    ! Add a minimum threshold to prevent degeneracy
    !mm_thresh = max(mm_thresh, 1.0d0)  ! Ensure mm_thresh is at least 1.0
    ! Initialize time_index = (1:m)/m
    do i = 1, m
        time_index(i) = dble(i) / dble(m)
    end do

    ! Build A matrix
    do i = 1, m
        A(i, 1) = 1.0d0
        do j = 1, nvar
            A(i, j+1) = x(i, j)
        end do
        do j = 1, nvar+1
            A(i, nvar+1+j) = A(i, j) * time_index(i)
        end do
    end do

    ! Initialize outputs
    theta_ll_est = 0.0d0
    beta_full_est = 0.0d0
    H_seq = 0
    full_active_recovery_count = 0
    if (threshold_scale_mode_int == 1) then
        threshold_scale = log(log(dble(m)))
    else if (threshold_scale_mode_int == 2) then
        threshold_scale = log(dble(m))
    else
        ierr = 5
        failed_eval = 1
        return
    end if
    if ((.not. threshold_lower_bound) .and. Mm_factor <= 0.0d0) then
        ierr = 5
        failed_eval = 1
        return
    end if

    ! ============================================
    ! EVA_T = 1: Standard simplex (no preprocessing)
    ! ============================================

    ! Calculate weights for t=1
    do i = 1, m
        if (abs(1.0d0/dble(m) - time_index(i)) <= h) then
            w(i) = 0.75d0 * (1.0d0 - ((1.0d0/dble(m) - time_index(i))/h)**2)
        else
            w(i) = 0.0d0
        end if
    end do

    ! Initialize cc vector
    cc = 0.0d0
    do i = 1, m
        cc(2*(1+nvar) + i) = tau * w(i)
        cc(2*(1+nvar) + m + i) = (1.0d0 - tau) * w(i)
    end do

    ! Build gammax matrix
    do i = 1, m
        do j = 1, 2*(nvar+1)
            gammax(i, j) = A(i, j)
        end do
    end do

    ! Initialize b
    do i = 1, m
        b(i) = y(i)
    end do
    b(m+1) = 0.0d0

    ! Adjust gammax and b based on sign of y - matching R logic
    do i = 1, m
        if (y(i) < 0.0d0) then
            do j = 1, 2*(nvar+1)
                gammax(i, j) = -gammax(i, j)
            end do
            b(i) = -b(i)
        end if
    end do

    ! Initialize IB - matching R: IB <- (y >= 0) * ((1:m) + 2 * (1 + nvar)) + (y < 0) * ((1:m) + 2 * (1 + nvar) + m)
    do i = 1, m
        if (y(i) >= 0.0d0) then
            IB(i) = i + 2*(1+nvar)
        else
            IB(i) = i + 2*(1+nvar) + m
        end if
    end do

    ! Add last row to gammax
    do j = 1, 2*(nvar+1)
        gammax(m+1, j) = 0.0d0
        do i = 1, m
            gammax(m+1, j) = gammax(m+1, j) - cc(IB(i)) * gammax(i, j)
        end do
    end do

    ! Initialize freevarrow
    freevarrow(1:m+1) = .false.
    freevarrow(m+1) = .true.

    ! Initialize r1 and r2
    do i = 1, 2*(nvar+1)
        r1(i) = i
        r2(i) = 0
    end do


    ! Simplex iterations for eva_t = 1
    iter = 0
    simplex_converged = .false.




    do while (iter < maxit)


! Compute reduced costs
        do i = 1, 2*(nvar+1)
            rr(1, i) = gammax(m+1, i)

            ! Calculate the index for weight lookup
            j = r1(i) - 2 - 2*nvar

            ! R logic: (r1 - 2 - 2 * nvar) * (r1 - 2 - 2 * nvar > 0) + (r1 - 2 - 2 * nvar <= 0)
            if (j > 0) then
                ! Use j when positive
                if (j <= m) then
                    rr(2, i) = w(j) - rr(1, i)
                else
                    rr(2, i) = w(1) - rr(1, i)  ! Bounds check
                end if
            else
                ! When j <= 0, the R expression evaluates to 1
                rr(2, i) = w(1) - rr(1, i)
            end if

            ! Apply the (r2 != 0) condition from R
            if (r2(i) == 0) then
                rr(2, i) = 0.0d0
            end if

            ! Handle the sign flip for r2 == 0
            if (r2(i) == 0) then
                rr(1, i) = -abs(rr(1, i))
            end if
        end do

        ! Check optimality
        rrl = minval(rr)
        if (rrl >= -tol) then
            simplex_converged = .true.
            exit
        end if

        ! Choose entering variable
        if (bland) then
            t = huge(1)
            t_rr = 0

            do i = 1, 2*(nvar+1)
                if (rr(1, i) < -tol) then
                    if (r1(i) < t) then
                        t = r1(i)
                        t_rr = i
                        tsep = 1
                    end if
                end if
            end do

            if (t_rr == 0) then
                do i = 1, 2*(nvar+1)
                    if (rr(2, i) < -tol) then
                        if (r2(i) < t) then
                            t = r2(i)
                            t_rr = i
                            tsep = 2
                        end if
                    end if
                end do
            end if
        else
            do j = 1, 2*(nvar+1)
                do i = 1, 2
                    if (abs(rr(i, j) - rrl) <= 1.0d-14) then
                        t_rr = j
                        tsep = i
                        if (tsep == 1) then
                            t = r1(t_rr)
                        else
                            t = r2(t_rr)
                        end if
                        goto 100
                    end if
                end do
            end do
            100 continue
        end if

        ! Choose leaving variable
        if (r2(t_rr) /= 0) then
            if (tsep == 1) then
                yy(1:m+1) = gammax(:, t_rr)
            else
                yy(1:m+1) = -gammax(:, t_rr)
            end if

            min_k = huge(1.0d0)
            k = 0

            do i = 1, m+1
                if (yy(i) > pivot_tol .and. .not. freevarrow(i)) then
                    k_vals(i) = b(i) / yy(i)
                    if (k_vals(i) < min_k) then
                        min_k = k_vals(i)
                        if (.not. bland) k = i
                    end if
                else
                    k_vals(i) = huge(1.0d0)
                end if
            end do

            if (bland .and. k == 0) then
                do i = 1, m+1
                    if (abs(k_vals(i) - min_k) < 1.0d-14) then
                        if (k == 0 .or. IB(i) < IB(k)) then
                            k = i
                        end if
                    end if
                end do
            end if

            if (k == 0 .or. min_k >= huge(1.0d0)) then
                call set_failure(2, 1)
                return
            end if

            if (tsep /= 1) then
                j = r1(t_rr) - 2 - 2*nvar
                if (j > 0 .and. j <= m) then
                    yy(m+1) = yy(m+1) + w(j)
                end if
            end if

        else
            yy(1:m+1) = gammax(:, t_rr)

            if (yy(m+1) < 0.0d0) then
                min_k = huge(1.0d0)
                k = 0

                do i = 1, m+1
                    if (yy(i) > pivot_tol .and. .not. freevarrow(i)) then
                        k_vals(i) = b(i) / yy(i)
                        if (k_vals(i) < min_k) then
                            min_k = k_vals(i)
                            if (.not. bland) k = i
                        end if
                    else
                        k_vals(i) = huge(1.0d0)
                    end if
                end do

                if (bland .and. k == 0) then
                    do i = 1, m+1
                        if (abs(k_vals(i) - min_k) < 1.0d-14) then
                            if (k == 0 .or. IB(i) < IB(k)) then
                                k = i
                            end if
                        end if
                    end do
                end if
            else
                min_k = huge(1.0d0)
                k = 0

                do i = 1, m+1
                    if (yy(i) < -pivot_tol .and. .not. freevarrow(i)) then
                        k_vals(i) = -b(i) / yy(i)
                        if (k_vals(i) < min_k) then
                            min_k = k_vals(i)
                            if (.not. bland) k = i
                        end if
                    else
                        k_vals(i) = huge(1.0d0)
                    end if
                end do

                if (bland .and. k == 0) then
                    do i = 1, m+1
                        if (abs(k_vals(i) - min_k) < 1.0d-14) then
                            if (k == 0 .or. IB(i) < IB(k)) then
                                k = i
                            end if
                        end if
                    end do
                end if
            end if

            if (k == 0 .or. min_k >= huge(1.0d0)) then
                call set_failure(2, 1)
                return
            end if

            freevarrow(k) = .true.
        end if

        ! Perform pivot
        ee(1:m+1) = yy(1:m+1) / yy(k)
        ee(k) = 1.0d0 - 1.0d0 / yy(k)

        if (IB(k) <= (m + 2*nvar + 2)) then
            gammax(:, t_rr) = 0.0d0
            gammax(k, t_rr) = 1.0d0
            r1(t_rr) = IB(k)
            r2(t_rr) = IB(k) + m
        else
            gammax(:, t_rr) = 0.0d0
            gammax(k, t_rr) = -1.0d0
            j = IB(k) - m - 2*nvar - 2
            if (j > 0 .and. j <= m) then
                gammax(m+1, t_rr) = w(j)
            else
                gammax(m+1, t_rr) = 0.0d0
            end if
            r1(t_rr) = IB(k) - m
            r2(t_rr) = IB(k)
        end if

        do j = 1, 2*(nvar+1)
            pivot_row(j) = gammax(k, j)
        end do

        do j = 1, 2*(nvar+1)
            do i = 1, m+1
                gammax(i, j) = gammax(i, j) - ee(i) * pivot_row(j)
            end do
        end do

        b_k_original = b(k)
        do i = 1, m+1
            b(i) = b(i) - ee(i) * b_k_original
        end do

        IB(k) = t

        iter = iter + 1
    end do

    if (.not. simplex_converged) then
        call set_failure(2, 1)
        return
    end if


    ! Extract solution for t=1
    u = 0.0d0
    v = 0.0d0
    estimate = 0.0d0




    do i = 1, m
        if (IB(i) >= 1 .and. IB(i) <= 2*(nvar+1)) then
            estimate(IB(i)) = b(i)
        else if (IB(i) >= 2*(nvar+1)+1 .and. IB(i) <= 2*(nvar+1)+m) then
            j = IB(i) - 2*(nvar+1)
            u(j) = b(i)
        else if (IB(i) > 2*(nvar+1)+m) then
            j = IB(i) - 2*(nvar+1) - m
            v(j) = b(i)
        end if
    end do

    do j = 1, nvar+1
        theta_ll_est(1, j) = estimate(j) + (1.0d0/dble(m)) * estimate(nvar+1+j)
    end do
    do j = 1, 2*(nvar+1)
        beta_full_est(1, j) = estimate(j)
    end do

    do i = 1, m
        r_prev(i) = u(i) - v(i)
    end do



    ! Extract H observations - only those that correspond to actual observations
    do i = 1, 2*(nvar+1)
        j = r1(i) - 2 - 2*nvar
        if (j > 0 .and. j <= m) then
            H_seq(1, i) = j
        else
            ! This r1 element doesn't correspond to an observation
            ! This can happen when a coefficient variable is in the basis
            H_seq(1, i) = 0  ! Mark as invalid
        end if
    end do
    do i = 1, 2*(nvar+1)
        if (H_seq(1, i) < 1 .or. H_seq(1, i) > m) then
            call set_failure(3, 1)
            return
        end if
    end do

    ! Save X(h)^{-1} for next iteration (stored in first 2*(nvar+1) rows of gammax)
    do i = 1, 2*(nvar+1)
        do j = 1, 2*(nvar+1)
            gammaxs_temp(i, j) = gammax(i, j)
        end do
    end do

    ! ============================================
    ! EVA_T = 2 to 4: Preprocessing iterations
    ! ============================================

    do eva_t = 2, m






        not_optimal = .true.
        not_new_sl_sh = .true.
        force_full_sample = .false.
        empty_pivot_count = 0









        ! Initialize counters for this time point
        !n_Hbar_pos = 0
        !n_Hbar_neg = 0

        ! Update weights for current t
        n_active = 0
        do i = 1, m
            if (abs(dble(eva_t)/dble(m) - time_index(i)) <= h) then
                w(i) = 0.75d0 * (1.0d0 - ((dble(eva_t)/dble(m) - time_index(i))/h)**2)
            else
                w(i) = 0.0d0
            end if
            active(i) = (w(i) > 0.0d0)
            if (active(i)) n_active = n_active + 1
        end do





        mmm_thresh = mm_thresh

        j = 0
        iter = 0  ! Initialize iteration counter here, outside preprocessing loop
        use_independent_init = .false.
        independent_trigger = 0

        preprocessing_attempts = 0
        do while (not_optimal)
            unbounded_detected = .false.  ! Initialize for each attempt
            simplex_converged = .false.
            no_pivot_flag = .false.
            iter_attempt = 0

            preprocessing_attempts = preprocessing_attempts + 1






            total_preprocessing_loops = total_preprocessing_loops + 1
            ! Add safety valve to prevent infinite preprocessing
            if (preprocessing_attempts > max(10, 2*(nvar+1))) then
                force_full_sample = .true.
                not_new_sl_sh = .true.
            end if
            ! Get residuals from previous time
            do i = 1, m
                r(i) = r_prev(i)
            end do

            if (.not. use_independent_init) then
                call validate_previous_H(H_seq(eva_t-1, :), previous_H_valid)
                if (.not. previous_H_valid) then
                    use_independent_init = .true.
                    independent_trigger = 1
                end if
            end if





            ! Calculate threshold and partition observations. sl/sh are only
            ! recomputed when requested; all dependent masks/counts are rebuilt
            ! below on every attempt.
            if (force_full_sample) then
                sl = .false.
                sh = .false.
            else if (not_new_sl_sh) then
                ! Create array of absolute residuals
                do i = 1, m
                    abs_r(i) = abs(r(i))
                end do
                ! Fix 4: Scale threshold to residual magnitude (KEEP THIS!)
                residual_scale = median_value(abs_r, m)
                if (threshold_lower_bound) then
                    M_threshold = max(Mm_factor * mmm_thresh * threshold_scale, 0.1d0 * residual_scale)
                else
                    M_threshold = Mm_factor * mmm_thresh * threshold_scale
                end if

                ! NEW: Ensure minimum subsample size
                if (min_subsample_size_in >= 0) then
                    min_subsample_size = min_subsample_size_in
                else
                    min_subsample_size = max(5 * (2 * (nvar + 1)), ceiling(0.2d0 * dble(m)))
                end if

                target_min = min(min_subsample_size, n_active)

                ! Count active potential observations in S
                n_potential_S = 0
                min_k = 0.0d0
                do i = 1, m
                    if (active(i)) then
                        if (abs(r(i)) > min_k) min_k = abs(r(i))
                    end if
                    if (active(i) .and. abs(r(i)) <= M_threshold) then
                        n_potential_S = n_potential_S + 1
                    end if
                end do

                ! If too few active observations would remain, increase M.
                do while (n_potential_S < target_min .and. M_threshold < min_k)
                    M_threshold = M_threshold * 1.5d0
                    n_potential_S = 0
                    do i = 1, m
                        if (active(i) .and. abs(r(i)) <= M_threshold) then
                            n_potential_S = n_potential_S + 1
                        end if
                    end do
                end do


                do i = 1, m
                    sl(i) = active(i) .and. r(i) < -M_threshold
                    sh(i) = active(i) .and. r(i) > M_threshold
                end do




            end if






            ! Always rebuild the retained-subsample mask from active current sl/sh.
            if (force_full_sample) then
                sl = .false.
                sh = .false.
                not_jl_or_jh = active
            else
                do i = 1, m
                    not_jl_or_jh(i) = active(i) .and. (.not. (sl(i) .or. sh(i)))
                end do
            end if

            ! Always force previous H observations into the retained subsample.
            ! Zero-weight previous-H rows are basis padding, not active screened rows.
            if (.not. use_independent_init) then
                do k = 1, 2*(nvar+1)
                    idx = H_seq(eva_t-1, k)
                    sl(idx) = .false.
                    sh(idx) = .false.
                    not_jl_or_jh(idx) = .true.
                end do
            end if

            ! Recompute all counts and retained indices after force-H.
            n_sl = 0
            n_sh = 0
            sum_w_sl = 0.0d0
            sum_w_sh = 0.0d0
            ms_org = 0
            do i = 1, m
                if (not_jl_or_jh(i)) then
                    ms_org = ms_org + 1
                    idx_not_jl_or_jh(ms_org) = i
                end if
                if (sl(i)) then
                    n_sl = n_sl + 1
                    sum_w_sl = sum_w_sl + w(i)
                end if
                if (sh(i)) then
                    n_sh = n_sh + 1
                    sum_w_sh = sum_w_sh + w(i)
                end if
            end do
            has_sl_agg = (sum_w_sl > 0.0d0)
            has_sh_agg = (sum_w_sh > 0.0d0)
            ms = ms_org

            if (ms_org < 2*(nvar+1)) then
                if (force_full_sample) then
                    call set_failure(5, eva_t)
                    return
                end if
                mmm_thresh = 2.0d0 * mmm_thresh
                not_new_sl_sh = .true.
                if (preprocessing_attempts >= max_empty_pivot_retries) then
                    force_full_sample = .true.
                end if
                cycle
            end if





            ms_org = ms  ! Store the original subsample size before adding aggregated obs


            ! Extract subsample weights
            do i = 1, ms
                ws(i) = w(idx_not_jl_or_jh(i))
            end do


            ! Initialize for eva_t = 2
            if (eva_t == 2 .or. use_independent_init) then
                do i = 1, ms
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(i, j) = A(idx_not_jl_or_jh(i), j)
                    end do
                    bs_temp(i) = y(idx_not_jl_or_jh(i))
                end do




            end if

            ! Add aggregated observations if JL is not empty
            if (has_sl_agg) then
                glob_wx = 0.0d0
                glob_wy = 0.0d0

                do i = 1, m
                    if (sl(i)) then
                        wsl(i) = w(i)
                        do j = 1, 2*(nvar+1)
                            glob_wx(j) = glob_wx(j) + A(i, j) * wsl(i)
                        end do
                        glob_wy = glob_wy + y(i) * wsl(i)
                    end if
                end do





                do j = 1, 2*(nvar+1)
                    gammaxs_temp(m + 1, j) = glob_wx(j)
                end do
                bs_temp(m + 1) = glob_wy
                ms = ms + 1
                ws(ms) = 1.0d0




            end if

            ! Add aggregated observations if JH is not empty
            if (has_sh_agg) then
                ghib_wx = 0.0d0
                ghib_wy = 0.0d0

                do i = 1, m
                    if (sh(i)) then
                        wsh(i) = w(i)
                        do j = 1, 2*(nvar+1)
                            ghib_wx(j) = ghib_wx(j) + A(i, j) * wsh(i)
                        end do
                        ghib_wy = ghib_wy + y(i) * wsh(i)
                    end if
                end do



                do j = 1, 2*(nvar+1)
                    gammaxs_temp(m + 2, j) = ghib_wx(j)
                end do
                bs_temp(m + 2) = ghib_wy
                ms = ms + 1
                ws(ms) = 1.0d0



            end if



            ! Find positive and negative residuals in subsample
            n_idpos = 0
            n_idneg = 0
            do i = 1, ms_org
                if (r(idx_not_jl_or_jh(i)) > 0.0d0) then
                    n_idpos = n_idpos + 1
                    idpos(n_idpos) = i
                else if (r(idx_not_jl_or_jh(i)) < 0.0d0) then
                    n_idneg = n_idneg + 1
                    idneg(n_idneg) = i
                end if
            end do

            ! Map H indices to subsample
            do i = 1, 2*(nvar+1)

                H_indices(i) = 0
                if (H_seq(eva_t-1, i) > 0 .and. H_seq(eva_t-1, i) <= m) then
                    ! Valid H observation - try to find it in subsample
                    do j = 1, ms_org
                        if (H_seq(eva_t-1, i) == idx_not_jl_or_jh(j)) then
                            H_indices(i) = j
                            exit
                        end if
                    end do
                end if
            end do






            ! Check if any H_indices are missing (important for degenerate cases)
            k = 0  ! Count valid H observations
            do i = 1, 2*(nvar+1)
                if (H_seq(eva_t-1, i) > 0 .and. H_seq(eva_t-1, i) <= m) then
                    if (H_indices(i) == 0) then
                    else
                        k = k + 1
                    end if
                end if
            end do

            ! Check if we have enough valid H observations
            if ((.not. use_independent_init) .and. k < 2*(nvar+1)) then
                use_independent_init = .true.
                independent_trigger = 2
            end if

            if (.not. use_independent_init) then
                call matrix_inverse_2p(A(H_seq(eva_t-1, :), :), xhinv, 2*(nvar+1), inv_info)
                if (inv_info /= 0) then
                    use_independent_init = .true.
                    independent_trigger = 3
                end if
            end if

            if (use_independent_init) then
                do i = 1, ms_org
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(i, j) = A(idx_not_jl_or_jh(i), j)
                    end do
                    bs_temp(i) = y(idx_not_jl_or_jh(i))
                end do
                do i = 1, 2*(nvar+1)
                    beta_prev(i) = beta_full_est(eva_t-1, i)
                end do
                call run_fresh_initialization(beta_prev, ms_org, ms, shifted_init_success)
                if (.not. shifted_init_success) then
                    full_active_recovery_count = full_active_recovery_count + 1
                    call run_full_active_recovery(full_recovery_success)
                    if (.not. full_recovery_success) then
                        call set_failure(2, eva_t)
                        return
                    end if
                end if
                simplex_converged = .true.
                no_pivot_flag = .false.
                goto 24430
            end if

            ! CRITICAL: Reset freevarrow for EVERY preprocessing attempt
            ! This must be done before building the tableau
            freevarrow(1:ms+1) = .false.
            do i = 1, 2*(nvar+1)
                freevarrow(i) = .true.
            end do
            freevarrow(ms+1) = .true.  ! objective row


            ! Remove H from idpos and idneg
            do i = 1, 2*(nvar+1)
                if (H_indices(i) > 0) then
                    do j = 1, n_idpos
                        if (idpos(j) == H_indices(i)) then
                            idpos(j) = 0
                        end if
                    end do
                    do j = 1, n_idneg
                        if (idneg(j) == H_indices(i)) then
                            idneg(j) = 0
                        end if
                    end do
                end if
            end do





            ! Compact idpos and idneg
            k = 0
            do i = 1, n_idpos
                if (idpos(i) /= 0) then
                    k = k + 1
                    idpos(k) = idpos(i)
                end if
            end do
            n_idpos = k

            k = 0
            do i = 1, n_idneg
                if (idneg(i) /= 0) then
                    k = k + 1
                    idneg(k) = idneg(i)
                end if
            end do
            n_idneg = k


            ! Set up u_in_IBs and v_in_IBs
            do i = 1, n_idpos
                u_in_IBs(i) = idpos(i) + 2*(nvar + 1)
            end do
            do i = 1, n_idneg
                v_in_IBs(i) = idneg(i) + 2*(nvar + 1) + ms
            end do

            ! Set up r1 and r2 for subsample
            do i = 1, 2*(nvar+1)
                r1(i) = H_indices(i) + 2*(nvar + 1)
                r2(i) = r1(i) + ms
            end do






            ! Now handle the initialization based on eva_t
            if (eva_t == 2) then
                ! Initialize IBs and freevarrow based on four cases
                if (has_sl_agg .and. has_sh_agg) then
                    ! Case 1: Both JL and JH non-empty
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + 2*ms - 1
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 2) = 2*(nvar+1) + ms

                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms - 1
                    Hbar_indices(n_idpos + n_idneg + 2) = ms

                    do i = 1, ms-2
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms-1, j) = gammaxs_temp(m+1, j)
                        gammaxs(ms, j) = gammaxs_temp(m+2, j)
                    end do
                    bs(ms-1) = bs_temp(m+1)
                    bs(ms) = bs_temp(m+2)

                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
                    lambda(n_idpos + n_idneg + 2) = tau


                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms-1) = .true.
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.





                else if (has_sl_agg) then
                    ! Case 2: Only JL non-empty
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + 2*ms

                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms

                    do i = 1, ms-1
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms, j) = gammaxs_temp(m+1, j)
                    end do
                    bs(ms) = bs_temp(m+1)

                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau

                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.

                else if (has_sh_agg) then
                    ! Case 3: Only JH non-empty
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + ms

                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms

                    do i = 1, ms-1
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms, j) = gammaxs_temp(m+2, j)
                    end do
                    bs(ms) = bs_temp(m+2)

                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = tau

                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.




                else
                    ! Case 4: Both JL and JH empty
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do

                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do

                    do i = 1, ms
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do

                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do

                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms+1) = .true.
                end if

                ! Compute initial tableau for t=2
                call matrix_inverse_2p(gammaxs(H_indices, :), xhinv, 2*(nvar+1), inv_info)

                if (inv_info /= 0) then
                    if (force_full_sample) then
                        call set_failure(4, eva_t)
                        return
                    end if
                    empty_pivot_count = empty_pivot_count + 1
                    mmm_thresh = 2.0d0 * mmm_thresh
                    not_new_sl_sh = .true.
                    if (empty_pivot_count >= max_empty_pivot_retries) then
                        force_full_sample = .true.
                    end if
                    cycle
                end if

                ! Compute Pxhbarxhinv = P @ gammaxs[Hbar, :] @ xhinv
                ! For positive residuals (u variables), P has +1
                do i = 1, n_idpos
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(i, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(i, j) = Pxhbarxhinv(i, j) + gammaxs(Hbar_indices(i), k) * xhinv(k, j)
                        end do
                    end do
                end do

                ! For negative residuals (v variables), P has -1
                do i = 1, n_idneg
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + i, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + i, j) = Pxhbarxhinv(n_idpos + i, j) - &
                                                         gammaxs(Hbar_indices(n_idpos + i), k) * xhinv(k, j)
                        end do
                    end do
                end do

                ! Handle JL and JH based on case
                if (has_sl_agg .and. has_sh_agg) then
                    ! For v_L (JL aggregated), P has -1
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + n_idneg + 1, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + n_idneg + 1, j) = Pxhbarxhinv(n_idpos + n_idneg + 1, j) - &
                                                                   gammaxs(ms-1, k) * xhinv(k, j)
                        end do
                    end do

                    ! For u_H (JH aggregated), P has +1
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + n_idneg + 2, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + n_idneg + 2, j) = Pxhbarxhinv(n_idpos + n_idneg + 2, j) + &
                                                                   gammaxs(ms, k) * xhinv(k, j)
                        end do
                    end do
                else if (has_sl_agg) then
                    ! Only v_L
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + n_idneg + 1, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + n_idneg + 1, j) = Pxhbarxhinv(n_idpos + n_idneg + 1, j) - &
                                                                   gammaxs(ms, k) * xhinv(k, j)
                        end do
                    end do
                else if (has_sh_agg) then
                    ! Only u_H
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(n_idpos + n_idneg + 1, j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            Pxhbarxhinv(n_idpos + n_idneg + 1, j) = Pxhbarxhinv(n_idpos + n_idneg + 1, j) + &
                                                                   gammaxs(ms, k) * xhinv(k, j)
                        end do
                    end do
                end if

                ! Build the tableau
                ! First rows are xhinv
                do i = 1, 2*(nvar+1)
                    do j = 1, 2*(nvar+1)
                        gammaxs(i, j) = xhinv(i, j)
                    end do
                end do




                ! Next rows are -Pxhbarxhinv
                k = n_idpos + n_idneg
                if (has_sl_agg) k = k + 1
                if (has_sh_agg) k = k + 1

                do i = 1, k
                    do j = 1, 2*(nvar+1)
                        gammaxs(2*(nvar+1) + i, j) = -Pxhbarxhinv(i, j)
                    end do
                end do

                ! Last row is the objective function row
                do j = 1, 2*(nvar+1)
                    if (H_indices(j) < 1 .or. H_indices(j) > ms) then
                        gammaxs(ms + 1, j) = 0.0d0  ! Safety default
                    else
                        gammaxs(ms + 1, j) = tau * ws(H_indices(j))
                    end if
                    do i = 1, k
                        gammaxs(ms + 1, j) = gammaxs(ms + 1, j) + lambda(i) * Pxhbarxhinv(i, j)
                    end do
                end do

                ! Build bs vector
                ! First part: xhinv @ bs[H]


                do i = 1, 2*(nvar+1)
                    bs(i) = 0.0d0
                    do j = 1, 2*(nvar+1)
                        bs(i) = bs(i) + xhinv(i, j) * bs_temp(H_indices(j))
                    end do
                end do

                ! Second part: -Pxhbarxhinv @ bs[H] + P @ bs[Hbar]
                ! For positive residuals
                do i = 1, n_idpos
                    bs(2*(nvar+1) + i) = bs_temp(Hbar_indices(i))
                    do j = 1, 2*(nvar+1)
                        bs(2*(nvar+1) + i) = bs(2*(nvar+1) + i) - Pxhbarxhinv(i, j) * bs_temp(H_indices(j))
                    end do
                end do

                ! For negative residuals
                do i = 1, n_idneg
                    bs(2*(nvar+1) + n_idpos + i) = -bs_temp(Hbar_indices(n_idpos + i))
                    do j = 1, 2*(nvar+1)
                        bs(2*(nvar+1) + n_idpos + i) = bs(2*(nvar+1) + n_idpos + i) - &
                                                      Pxhbarxhinv(n_idpos + i, j) * bs_temp(H_indices(j))
                    end do
                end do

                ! Handle JL and JH
                if (has_sl_agg .and. has_sh_agg) then
                    bs(ms-1) = -bs_temp(m+1)
                    do j = 1, 2*(nvar+1)
                        bs(ms-1) = bs(ms-1) - Pxhbarxhinv(n_idpos + n_idneg + 1, j) * bs_temp(H_indices(j))
                    end do

                    bs(ms) = bs_temp(m+2)
                    do j = 1, 2*(nvar+1)
                        bs(ms) = bs(ms) - Pxhbarxhinv(n_idpos + n_idneg + 2, j) * bs_temp(H_indices(j))
                    end do
                else if (has_sl_agg) then
                    bs(ms) = -bs_temp(m+1)
                    do j = 1, 2*(nvar+1)
                        bs(ms) = bs(ms) - Pxhbarxhinv(n_idpos + n_idneg + 1, j) * bs_temp(H_indices(j))
                    end do
                else if (has_sh_agg) then
                    bs(ms) = bs_temp(m+2)
                    do j = 1, 2*(nvar+1)
                        bs(ms) = bs(ms) - Pxhbarxhinv(n_idpos + n_idneg + 1, j) * bs_temp(H_indices(j))
                    end do
                end if

                bs(ms+1) = 0.0d0




            else
                ! eva_t >= 3: Reuse computations from previous iteration





                ! xhinv is already stored in gammaxs_temp from previous iteration



                do i = 1, 2*(nvar+1)
                    do j = 1, 2*(nvar+1)
                        xhinv(i, j) = gammaxs_temp(i, j)
                    end do
                end do




                ! Map current idpos to original indices
                do i = 1, n_idpos
                    idx_Hbar_pos2(i) = idx_not_jl_or_jh(idpos(i))
                end do

                ! Check which positive residuals can be reused
                n_matched_pos = 0
                do i = 1, n_idpos
                    matched_rows_pos(i) = .false.
                    do j = 1, n_Hbar_pos
                        if (idx_Hbar_pos2(i) == idx_Hbar_pos(j)) then
                            matched_rows_pos(i) = .true.
                            n_matched_pos = n_matched_pos + 1
                            rows_pos(n_matched_pos) = j
                            exit
                        end if
                    end do
                end do





                ! Initialize positive residual rows in gammaxs_temp
                ! First, copy the rows that can be reused
                k = 0
                do i = 1, n_idpos
                    if (matched_rows_pos(i)) then
                        k = k + 1
                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(2*(nvar+1) + i, j) = gammaxs_pos(rows_pos(k), j)
                        end do
                        bs_temp(2*(nvar+1) + i) = bs_pos(rows_pos(k))
                    end if
                end do

                ! For new positive residuals, compute fresh
                do i = 1, n_idpos
                    if (.not. matched_rows_pos(i)) then
                        ! Compute -P @ X(hbar) @ X(h)^{-1}
                        do j = 1, 2*(nvar+1)
                            temp_vec(j) = 0.0d0
                            do k = 1, 2*(nvar+1)
                                temp_vec(j) = temp_vec(j) + A(idx_Hbar_pos2(i), k) * xhinv(k, j)
                            end do
                        end do

                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(2*(nvar+1) + i, j) = -temp_vec(j)
                        end do

                        ! Compute -P @ X(hbar) @ X(h)^{-1} @ y(h) + P @ y(hbar)
                        bs_temp(2*(nvar+1) + i) = y(idx_Hbar_pos2(i))
                        do j = 1, 2*(nvar+1)
                            bs_temp(2*(nvar+1) + i) = bs_temp(2*(nvar+1) + i) - &
                                                    temp_vec(j) * y(idx_not_jl_or_jh(H_indices(j)))
                        end do
                    end if
                end do

                ! Similarly for negative residuals
                do i = 1, n_idneg
                    idx_Hbar_neg2(i) = idx_not_jl_or_jh(idneg(i))
                end do

                n_matched_neg = 0
                do i = 1, n_idneg
                    matched_rows_neg(i) = .false.
                    do j = 1, n_Hbar_neg
                        if (idx_Hbar_neg2(i) == idx_Hbar_neg(j)) then
                            matched_rows_neg(i) = .true.
                            n_matched_neg = n_matched_neg + 1
                            rows_neg(n_matched_neg) = j
                            exit
                        end if
                    end do
                end do










                ! Copy reusable negative residual rows
                k = 0
                do i = 1, n_idneg
                    if (matched_rows_neg(i)) then
                        k = k + 1



                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(2*(nvar+1) + n_idpos + i, j) = gammaxs_neg(rows_neg(k), j)
                        end do
                        bs_temp(2*(nvar+1) + n_idpos + i) = bs_neg(rows_neg(k))
                    end if
                end do




                ! Compute fresh for new negative residuals
                do i = 1, n_idneg
                    if (.not. matched_rows_neg(i)) then
                        ! Compute P @ X(hbar) @ X(h)^{-1} (note: P is -1 for v variables)
                        do j = 1, 2*(nvar+1)
                            temp_vec(j) = 0.0d0
                            do k = 1, 2*(nvar+1)
                                temp_vec(j) = temp_vec(j) - A(idx_Hbar_neg2(i), k) * xhinv(k, j)
                            end do
                        end do

                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(2*(nvar+1) + n_idpos + i, j) = -temp_vec(j)
                        end do

                        bs_temp(2*(nvar+1) + n_idpos + i) = -y(idx_Hbar_neg2(i))
                        do j = 1, 2*(nvar+1)
                            bs_temp(2*(nvar+1) + n_idpos + i) = bs_temp(2*(nvar+1) + n_idpos + i) - &
                                                              temp_vec(j) * y(idx_not_jl_or_jh(H_indices(j)))
                        end do
                    end if
                end do



                ! Handle JL and JH aggregated rows
                if (has_sl_agg .and. has_sh_agg) then



                    ! For JL: compute X_L^T @ X(h)^{-1}
                    do j = 1, 2*(nvar+1)
                        temp_vec(j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            temp_vec(j) = temp_vec(j) + gammaxs_temp(m+1, k) * xhinv(k, j)
                        end do

                    end do

                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(m+1, j) = temp_vec(j)
                    end do



                    ! Compute corresponding bs entry
                    temp_sum = 0.0d0
                    do j = 1, 2*(nvar+1)
                        temp_sum = temp_sum + gammaxs_temp(m+1, j) * y(idx_not_jl_or_jh(H_indices(j)))
                    end do
                    bs_temp(m+1) = temp_sum - bs_temp(m+1)



                    ! For JH: compute -X_H^T @ X(h)^{-1}
                    do j = 1, 2*(nvar+1)
                        temp_vec(j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            temp_vec(j) = temp_vec(j) - gammaxs_temp(m+2, k) * xhinv(k, j)
                        end do
                        !gammaxs_temp(m+2, j) = temp_vec(j)
                    end do

                    ! Now copy the result back AFTER computing all elements
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(m+2, j) = temp_vec(j)
                    end do

                    temp_sum = 0.0d0
                    do j = 1, 2*(nvar+1)
                        temp_sum = temp_sum + gammaxs_temp(m+2, j) * y(idx_not_jl_or_jh(H_indices(j)))
                    end do
                    bs_temp(m+2) = temp_sum + bs_temp(m+2)

                    ! Copy to gammaxs
                    do i = 1, ms-2
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms-1, j) = gammaxs_temp(m+1, j)
                        gammaxs(ms, j) = gammaxs_temp(m+2, j)
                    end do
                    bs(ms-1) = bs_temp(m+1)
                    bs(ms) = bs_temp(m+2)

                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
                    lambda(n_idpos + n_idneg + 2) = tau

                else if (has_sl_agg) then
                    ! Similar for only JL case
                    do j = 1, 2*(nvar+1)
                        temp_vec(j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            temp_vec(j) = temp_vec(j) + gammaxs_temp(m+1, k) * xhinv(k, j)
                        end do

                    end do

                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(m+1, j) = temp_vec(j)
                    end do






                    temp_sum = 0.0d0
                    do j = 1, 2*(nvar+1)
                        temp_sum = temp_sum + gammaxs_temp(m+1, j) * y(idx_not_jl_or_jh(H_indices(j)))
                    end do
                    bs_temp(m+1) = temp_sum - bs_temp(m+1)

                    do i = 1, ms-1
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms, j) = gammaxs_temp(m+1, j)
                    end do
                    bs(ms) = bs_temp(m+1)




                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau



                else if (has_sh_agg) then





                    ! Similar for only JH case
                    do j = 1, 2*(nvar+1)
                        temp_vec(j) = 0.0d0
                        do k = 1, 2*(nvar+1)
                            temp_vec(j) = temp_vec(j) - gammaxs_temp(m+2, k) * xhinv(k, j)
                        end do
                        !gammaxs_temp(m+2, j) = temp_vec(j)
                    end do

                    ! Now copy the result back AFTER computing all elements
                    do j = 1, 2*(nvar+1)
                        gammaxs_temp(m+2, j) = temp_vec(j)
                    end do





                    temp_sum = 0.0d0
                    do j = 1, 2*(nvar+1)
                        temp_sum = temp_sum + gammaxs_temp(m+2, j) * y(idx_not_jl_or_jh(H_indices(j)))
                    end do
                    bs_temp(m+2) = temp_sum + bs_temp(m+2)

                    do i = 1, ms-1
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do
                    do j = 1, 2*(nvar+1)
                        gammaxs(ms, j) = gammaxs_temp(m+2, j)
                    end do
                    bs(ms) = bs_temp(m+2)

                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                    lambda(n_idpos + n_idneg + 1) = tau

                else
                    ! Both empty
                    do i = 1, ms
                        do j = 1, 2*(nvar+1)
                            gammaxs(i, j) = gammaxs_temp(i, j)
                        end do
                        bs(i) = bs_temp(i)
                    end do

                    do i = 1, n_idpos
                        lambda(i) = tau * ws(idpos(i))
                    end do
                    do i = 1, n_idneg
                        lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                    end do
                end if

                ! Set up IBs and freevarrow (same logic as eva_t = 2)
                if (has_sl_agg .and. has_sh_agg) then
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + 2*ms - 1
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 2) = 2*(nvar+1) + ms

                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms - 1
                    Hbar_indices(n_idpos + n_idneg + 2) = ms

                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms-1) = .true.
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.





                else if (has_sl_agg) then
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + 2*ms

                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms

                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.

                else if (has_sh_agg) then
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do
                    IBs(2*(nvar+1) + n_idpos + n_idneg + 1) = 2*(nvar+1) + ms

                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do
                    Hbar_indices(n_idpos + n_idneg + 1) = ms

                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms) = .true.
                    freevarrow(ms+1) = .true.

                else
                    do i = 1, 2*(nvar+1)
                        IBs(i) = i
                    end do
                    do i = 1, n_idpos
                        IBs(2*(nvar+1) + i) = u_in_IBs(i)
                    end do
                    do i = 1, n_idneg
                        IBs(2*(nvar+1) + n_idpos + i) = v_in_IBs(i)
                    end do

                    do i = 1, n_idpos
                        Hbar_indices(i) = idpos(i)
                    end do
                    do i = 1, n_idneg
                        Hbar_indices(n_idpos + i) = idneg(i)
                    end do

                    freevarrow(1:ms+1) = .false.
                    do i = 1, 2*(nvar+1)
                        freevarrow(i) = .true.
                    end do
                    freevarrow(ms+1) = .true.
                end if

                ! Extract Pxhbarxhinv from the stored -Pxhbarxhinv in gammaxs
                ! R: Pxhbarxhinv <- - gammaxs[(2*(nvar+1)+1):ms,]
                n_hbar_rows = ms - 2*(nvar+1)  ! This is the number of Hbar rows




                do i = 1, n_hbar_rows  ! Use the descriptive name
                    do j = 1, 2*(nvar+1)
                        Pxhbarxhinv(i, j) = -gammaxs(2*(nvar+1) + i, j)
                    end do
                end do







                ! Build the last row
                do j = 1, 2*(nvar+1)

                    ! Add bounds checking for H_indices
                    if (H_indices(j) > 0 .and. H_indices(j) <= ms_org) then
                        gammaxs(ms + 1, j) = tau * ws(H_indices(j))
                    else
                        ! This should not happen if H observations are correctly forced into subsample
                        call set_failure(3, eva_t)
                        return
                    end if






                    do i = 1, n_hbar_rows  ! Use the descriptive name
                        gammaxs(ms + 1, j) = gammaxs(ms + 1, j) + lambda(i) * Pxhbarxhinv(i, j)
                    end do
                end do








                bs(ms+1) = 0.0d0



            end if  ! End of eva_t == 2 vs eva_t >= 3




            ! Invariant checks before simplex
            if (ms < 2*(nvar+1)) then
                if (force_full_sample) then
                    call set_failure(5, eva_t)
                    return
                end if
                empty_pivot_count = empty_pivot_count + 1
                mmm_thresh = 2.0d0 * mmm_thresh
                not_new_sl_sh = .true.
                if (empty_pivot_count >= max_empty_pivot_retries) then
                    force_full_sample = .true.
                end if
                cycle
            end if

            ! Check that all beta columns are basic
            if (.not. all(freevarrow(1:2*(nvar+1)))) then
                do i = 1, 2*(nvar+1)
                    if (.not. freevarrow(i)) then
                    end if
                end do
            end if





            ! Simplex iterations for the reduced problem
            !iter = 0



            do while (iter < maxit)
                total_preprocessing_loops = total_preprocessing_loops + 1
                ! Add a safety check to prevent infinite preprocessing:
                if (preprocessing_attempts > 100) then
                    call set_failure(2, eva_t)
                    return
                end if
                ! Step 2: Compute reduced costs - CRITICAL SECTION
                ! This must match R: rr[2, ] <- (ws[r1 - 2 - 2 * nvar] - rr[1, ])
                do i = 1, 2*(nvar+1)
                      rr(1, i) = gammaxs(ms + 1, i)



                      ! Get the subsample index
                      j = r1(i) - 2 - 2*nvar

                      ! Use subsample weights ws, matching R's behavior
                      if (j > 0 .and. j <= ms) then
                          rr(2, i) = ws(j) - rr(1, i)
                      else if (j <= 0) then
                          ! R would use ws[1] when index is non-positive
                          rr(2, i) = ws(1) - rr(1, i)
                      else
                        ! j > ms: use ws(1), matching the reference behavior
                          rr(2, i) = ws(1) - rr(1, i)
                      end if
                end do




                if (eva_t == 1) then
                    do i = 1, 2*(nvar+1)
                        if (r2(i) == 0) then
                            rr(2, i) = 0.0d0
                            rr(1, i) = -abs(rr(1, i))
                        end if
                    end do
                end if





                ! Check optimality
                rrl = minval(rr)


                ! Verify the mask is working correctly
                if (rrl >= 0.0d0) then
                    do i = 1, 2*(nvar+1)
                        if (r2(i) == 0 .and. rr(2, i) /= 0.0d0) then
                            call set_failure(5, eva_t)
                            return
                        end if
                    end do
                end if

                ! Find location of minimum
                call minloc2d(rr, rrl, tsep, t_rr)




                ! ADD SAFETY VALVE HERE:
                if (iter > 50000) then



                    no_pivot_flag = .true.
                    exit  ! Break out of the simplex loop
                end if







                if (rrl >= -tol) then




                    simplex_converged = .true.  ! Mark as converged
                    exit
                end if




                ! Step 3: Choose entering variable
                if (bland) then
                    t = huge(1)
                    t_rr = 0

                    ! Check first row of reduced costs
                    do i = 1, 2*(nvar+1)
                        if (rr(1, i) < -tol) then
                            if (r1(i) < t) then
                                t = r1(i)
                                t_rr = i
                                tsep = 1
                            end if
                        end if
                    end do

                    ! If nothing found in first row, check second row
                    if (t_rr == 0) then
                        do i = 1, 2*(nvar+1)
                            if (rr(2, i) < -tol) then
                                if (r2(i) < t) then
                                    t = r2(i)
                                    t_rr = i
                                    tsep = 2
                                end if
                            end if
                        end do
                    end if
                else
                    ! Steepest descent rule
                    do j = 1, 2*(nvar+1)
                        do i = 1, 2
                            if (abs(rr(i, j) - rrl) < 1.0d-14) then
                                t_rr = j
                                tsep = i
                                if (tsep == 1) then
                                    t = r1(t_rr)
                                else
                                    t = r2(t_rr)
                                end if
                                goto 200
                            end if
                        end do
                    end do
                    200 continue
                end if





                ! Step 4: Get pivot column
                if (tsep == 1) then
                    yy(1:ms+1) = gammaxs(1:ms+1, t_rr)
                else
                    yy(1:ms+1) = -gammaxs(1:ms+1, t_rr)
                end if




                ! Step 5: Choose leaving variable (ratio test)
                min_k = huge(1.0d0)
                k = 0


                ! First, find all valid ratios
                do i = 1, ms
                    if (yy(i) > pivot_tol .and. .not. freevarrow(i)) then
                        k_vals(i) = bs(i) / yy(i)
                        if (k_vals(i) < min_k) then  !k_vals(i) >= -tol .and.
                            min_k = k_vals(i)
                            if (.not. bland) then
                                k = i
                            end if
                        end if
                    else
                        k_vals(i) = huge(1.0d0)
                    end if
                end do





                ! Check if we found any valid pivots
                if (min_k >= huge(1.0d0)) then







                    ! No valid pivots found - problem is unbounded



                    no_pivot_flag = .true.
                    unbounded_detected = .true.  ! SET THE FLAG
                    ! Exit the simplex loop to restart preprocessing
                    exit
                end if




                ! For Bland rule, choose the one with smallest index
                if (bland .and. min_k < huge(1.0d0)) then
                    k = 0
                    do i = 1, ms
                        if (abs(k_vals(i) - min_k) < 1.0d-14) then
                            if (k == 0 .or. IBs(i) < IBs(k)) then
                                k = i
                            end if
                        end if
                    end do
                end if




                ! Check if we found a valid pivot
                if (k == 0 .or. min_k >= huge(1.0d0)) then
                    no_pivot_flag = .true.
                    unbounded_detected = .true.
                    exit
                end if

                ! Pivoting step 6': Adjust for entering v_i
                if (tsep /= 1) then
                    ! Use subsample weight exactly as R does
                    j = r1(t_rr) - 2 - 2*nvar
                    if (j >= 1) then
                        yy(ms + 1) = yy(ms + 1) + ws(j)
                    end if
                end if


                ! Check if we actually found a valid pivot
                if (k == 0) then
                    no_pivot_flag = .true.
                    unbounded_detected = .true.
                    exit
                end if


                ! Perform pivot only when the selected pivot is safely away from zero.
                if (abs(yy(k)) <= pivot_tol) then
                    no_pivot_flag = .true.
                    exit
                end if

                ee(1:ms+1) = yy(1:ms+1) / yy(k)
                ee(k) = 1.0d0 - 1.0d0 / yy(k)





                ! Update basis representation
                if (IBs(k) <= (ms + 2*nvar + 2)) then
                    gammaxs(:, t_rr) = 0.0d0
                    gammaxs(k, t_rr) = 1.0d0
                    r1(t_rr) = IBs(k)
                    r2(t_rr) = IBs(k) + ms
                else
                    gammaxs(:, t_rr) = 0.0d0
                    gammaxs(k, t_rr) = -1.0d0
                    j = IBs(k) - ms - 2*nvar - 2
                    if (j >= 1 .and. j <= ms) then
                        gammaxs(ms + 1, t_rr) = ws(j)
                    else
                        gammaxs(ms + 1, t_rr) = 0.0d0
                    end if
                    r1(t_rr) = IBs(k) - ms
                    r2(t_rr) = IBs(k)
                end if

                ! Update tableau
                do j = 1, 2*(nvar+1)
                    pivot_row(j) = gammaxs(k, j)
                end do





                do j = 1, 2*(nvar+1)
                    do i = 1, ms+1
                        gammaxs(i, j) = gammaxs(i, j) - ee(i) * pivot_row(j)
                    end do
                end do

                b_k_original = bs(k)
                do i = 1, ms+1
                    bs(i) = bs(i) - ee(i) * b_k_original
                end do

                IBs(k) = t
                iter = iter + 1
                iter_attempt = iter_attempt + 1



            end do




            ! ADD ITERATION TRACKING HERE:
            total_simplex_iterations = total_simplex_iterations + iter_attempt
            if (iter > max_iter_at_any_t) then
                max_iter_at_any_t = iter
            end if

            ! Log concerning patterns:
            if (iter > 10000) then
            end if
            if (no_pivot_flag .or. (.not. simplex_converged)) then
                if (force_full_sample) then
                    call set_failure(2, eva_t)
                    return
                end if
                empty_pivot_count = empty_pivot_count + 1
                mmm_thresh = 2.0d0 * mmm_thresh
                not_new_sl_sh = .true.
                if (empty_pivot_count >= max_empty_pivot_retries .or. iter >= maxit) then
                    force_full_sample = .true.
                end if
                cycle
            end if

24430       continue
            do i = 1, 2*(nvar+1)
                estimate(i) = bs(i)
            end do

            ! Extract the final H from the converged reduced simplex basis.
            do i = 1, 2*(nvar+1)
                j = r1(i) - 2 - 2*nvar
                if (j > 0 .and. j <= ms_org) then
                    H_indices(i) = idx_not_jl_or_jh(j)
                else
                    H_indices(i) = 0
                end if
            end do

            refit_used = .false.
            if (always_same_h_refit .and. tvcqr_H_basis_valid(H_indices, A, m, 2*(nvar+1))) then
                call same_h_refit_tvcqr(H_indices, A, y, m, 2*(nvar+1), res_tol, &
                                        estimate_refit, r_refit, same_h_ok)
                if (same_h_ok) then
                    do i = 1, 2*(nvar+1)
                        estimate(i) = estimate_refit(i)
                    end do
                    do i = 1, m
                        r(i) = r_refit(i)
                    end do
                    refit_used = .true.
                end if
            end if

            ! Check signs of residuals using FULL sample. If same-H refit
            ! failed or was not attempted, fall back to the tableau estimate.
            if (.not. refit_used) then
                do i = 1, m
                    r(i) = y(i)
                    do j = 1, 2*(nvar+1)
                        r(i) = r(i) - A(i, j) * estimate(j)
                    end do
                end do
            end if




            ! Count bad signs using the full residual from the current candidate.
            bad_signs = 0
            n_sure_signs = n_sl + n_sh  ! Total sure-sign observations

            do i = 1, m
                if ((r(i) <= 0.0d0) .and. sh(i)) then
                    bad_signs = bad_signs + 1
                end if
                if ((r(i) >= 0.0d0) .and. sl(i)) then
                    bad_signs = bad_signs + 1
                end if
            end do




            ! Handle bad signs
            if (bad_signs > 0) then
                call handle_current_bad_signs()

            else
                ! No bad signs - we've reached optimality
                accept_subsample = certify_tvcqr_candidate(H_indices, r, A, m, 2*(nvar+1), res_tol)
                if ((.not. accept_subsample) .and. (.not. always_same_h_refit)) then
                    if (tvcqr_H_basis_valid(H_indices, A, m, 2*(nvar+1))) then
                        call same_h_refit_tvcqr(H_indices, A, y, m, 2*(nvar+1), res_tol, &
                                                estimate_refit, r_refit, same_h_ok)
                        if (same_h_ok) then
                            do i = 1, 2*(nvar+1)
                                estimate(i) = estimate_refit(i)
                            end do
                            do i = 1, m
                                r(i) = r_refit(i)
                            end do
                            refit_used = .true.
                            bad_signs = 0
                            do i = 1, m
                                if ((r(i) <= 0.0d0) .and. sh(i)) then
                                    bad_signs = bad_signs + 1
                                end if
                                if ((r(i) >= 0.0d0) .and. sl(i)) then
                                    bad_signs = bad_signs + 1
                                end if
                            end do
                            if (bad_signs == 0) then
                                accept_subsample = certify_tvcqr_candidate(H_indices, r, A, m, &
                                                                           2*(nvar+1), res_tol)
                            end if
                        end if
                    end if
                end if
                if (bad_signs > 0) then
                    call handle_current_bad_signs()
                    cycle
                end if
                if (.not. accept_subsample) then
                    if (((.not. any(sl)) .and. (.not. any(sh))) .or. force_full_sample) then
                        call set_failure(1, eva_t)
                        return
                    end if
                    mmm_thresh = 2.0d0 * mmm_thresh
                    not_new_sl_sh = .true.
                    cycle
                end if
                not_optimal = .false.
                do i = 1, 2*(nvar+1)
                    bs(i) = estimate(i)
                end do
                do i = 1, 2*(nvar+1)
                    r(H_indices(i)) = 0.0d0
                end do






                ! Extract H for next iteration
                do i = 1, 2*(nvar+1)
                    if (H_indices(i) > 0) then
                        H_seq(eva_t, i) = H_indices(i)
                    else
                        ! This matches R's behavior - r1 should always give valid indices
                        ! If not, there's a bug in the simplex algorithm
                        call set_failure(3, eva_t)
                        return
                    end if




                end do






                ! If not the last time point, save information for reuse
                if (eva_t < m) then







                    ! Calculate ms_org for saving
                    if (has_sl_agg .and. has_sh_agg) then
                        ms_org = ms - 2
                    else if (has_sl_agg .or. has_sh_agg) then
                        ms_org = ms - 1
                    else
                        ms_org = ms
                    end if

                    do i = 1, 2*(nvar+1)
                        do j = 1, 2*(nvar+1)
                            gammaxs_temp(i, j) = gammaxs(i, j)
                        end do
                        bs_temp(i) = bs(i)
                    end do

                    ! Initialize counters
                    n_Hbar_pos = 0
                    n_Hbar_neg = 0

                    ! Only process Hbar information if there are non-aggregated observations
                    if (ms_org > 0) then






                        ! Process the IBs to extract Hbar information
                        do idx = 2*(nvar+1) + 1, ms_org
                            i = idx - 2*(nvar+1)! idx = 2*(nvar+1) + i
                            if (IBs(idx) > 2*(nvar+1)+ms) then
                                ! This is a v variable
                                n_Hbar_neg = n_Hbar_neg + 1






                                id_gammaxs_Hbar(i) = IBs(idx) - 2*(nvar+1) - ms
                                if (id_gammaxs_Hbar(i) < 1 .or. id_gammaxs_Hbar(i) > ms_org) then
                                    call set_failure(5, eva_t)
                                    return
                                end if
                                idx_Hbar_neg(n_Hbar_neg) = idx_not_jl_or_jh(id_gammaxs_Hbar(i))
                                if (idx_Hbar_neg(n_Hbar_neg) < 1 .or. idx_Hbar_neg(n_Hbar_neg) > m) then
                                    call set_failure(5, eva_t)
                                    return
                                end if
                                do j = 1, 2*(nvar+1)
                                    gammaxs_neg(n_Hbar_neg, j) = gammaxs(idx, j)
                                end do
                                bs_neg(n_Hbar_neg) = max(-r(idx_Hbar_neg(n_Hbar_neg)), 0.0d0)
                            else
                                ! This is a u variable
                                n_Hbar_pos = n_Hbar_pos + 1
                                id_gammaxs_Hbar(i) = IBs(idx) - 2*(nvar+1)
                                if (id_gammaxs_Hbar(i) < 1 .or. id_gammaxs_Hbar(i) > ms_org) then
                                    call set_failure(5, eva_t)
                                    return
                                end if
                                idx_Hbar_pos(n_Hbar_pos) = idx_not_jl_or_jh(id_gammaxs_Hbar(i))
                                if (idx_Hbar_pos(n_Hbar_pos) < 1 .or. idx_Hbar_pos(n_Hbar_pos) > m) then
                                    call set_failure(5, eva_t)
                                    return
                                end if
                                do j = 1, 2*(nvar+1)
                                    gammaxs_pos(n_Hbar_pos, j) = gammaxs(idx, j)
                                end do
                                bs_pos(n_Hbar_pos) = max(r(idx_Hbar_pos(n_Hbar_pos)), 0.0d0)
                            end if
                        end do




                    end if

                end if

            end if

        end do  ! End while not_optimal

        ! Store results for this time point


        ! Compute theta_ll_est for current eva_t
        do j = 1, nvar+1
            theta_ll_est(eva_t, j) = estimate(j) + (dble(eva_t)/dble(m)) * estimate(nvar+1+j)
        end do
        do j = 1, 2*(nvar+1)
            beta_full_est(eva_t, j) = estimate(j)
        end do







        do i = 1, m
            r_prev(i) = r(i)
        end do



    end do  ! End of eva_t loop


    deallocate(A)
    deallocate(gammax)
    deallocate(gammaxs_temp)
    deallocate(bs_temp)
    deallocate(gammaxs)
    deallocate(bs)

contains

    subroutine validate_previous_H(H_idx, valid)
        implicit none
        integer, intent(in) :: H_idx(2*(nvar+1))
        logical, intent(out) :: valid
        integer :: vi, vj

        valid = .true.
        do vi = 1, 2*(nvar+1)
            if (H_idx(vi) < 1 .or. H_idx(vi) > m) then
                valid = .false.
                return
            end if
            do vj = vi + 1, 2*(nvar+1)
                if (H_idx(vi) == H_idx(vj)) then
                    valid = .false.
                    return
                end if
            end do
        end do
    end subroutine validate_previous_H

    subroutine reorder_fresh_tableau(theta_offset, n_individual, n_total, success)
        implicit none
        double precision, intent(in) :: theta_offset(2*(nvar+1))
        integer, intent(in) :: n_individual, n_total
        logical, intent(out) :: success
        double precision :: gx_tmp(m+3, 2*(nvar+1)), bv_tmp(m+3)
        integer :: IB_tmp(m+3)
        logical :: fvr_tmp(m+3), used_row(m+2)
        integer :: ci, cj, src, dest, qdim

        success = .false.
        qdim = 2*(nvar+1)
        used_row = .false.

        do ci = 1, qdim
            src = 0
            do cj = 1, n_total
                if (IBs(cj) == ci) then
                    src = cj
                    exit
                end if
            end do
            if (src < 1 .or. src > n_individual) return
            used_row(src) = .true.
            do cj = 1, qdim
                gx_tmp(ci, cj) = gammaxs(src, cj)
            end do
            bv_tmp(ci) = bs(src) + theta_offset(ci)
            IB_tmp(ci) = IBs(src)
            fvr_tmp(ci) = freevarrow(src)
        end do

        dest = qdim
        do src = 1, n_individual
            if (.not. used_row(src)) then
                dest = dest + 1
                do cj = 1, qdim
                    gx_tmp(dest, cj) = gammaxs(src, cj)
                end do
                bv_tmp(dest) = bs(src)
                IB_tmp(dest) = IBs(src)
                fvr_tmp(dest) = freevarrow(src)
            end if
        end do
        do src = n_individual + 1, n_total
            dest = dest + 1
            do cj = 1, qdim
                gx_tmp(dest, cj) = gammaxs(src, cj)
            end do
            bv_tmp(dest) = bs(src)
            IB_tmp(dest) = IBs(src)
            fvr_tmp(dest) = freevarrow(src)
        end do
        if (dest /= n_total) return

        do cj = 1, qdim
            gx_tmp(n_total + 1, cj) = gammaxs(n_total + 1, cj)
        end do
        bv_tmp(n_total + 1) = bs(n_total + 1)
        IB_tmp(n_total + 1) = IBs(n_total + 1)
        fvr_tmp(n_total + 1) = freevarrow(n_total + 1)

        do ci = 1, n_total + 1
            do cj = 1, qdim
                gammaxs(ci, cj) = gx_tmp(ci, cj)
            end do
            bs(ci) = bv_tmp(ci)
            IBs(ci) = IB_tmp(ci)
            freevarrow(ci) = fvr_tmp(ci)
        end do
        success = .true.
    end subroutine reorder_fresh_tableau

    subroutine run_fresh_initialization(theta_offset, n_individual, n_total, success)
        implicit none
        double precision, intent(in) :: theta_offset(2*(nvar+1))
        integer, intent(in) :: n_individual, n_total
        logical, intent(out) :: success
        double precision :: shifted_rhs, delta_est(2*(nvar+1))
        integer :: si, sj, src, qdim, reduced_idx, iter_fresh, remaining
        logical :: no_pivot_fresh, converged_fresh, reorder_ok

        success = .false.
        qdim = 2*(nvar+1)
        if (n_individual < qdim .or. n_total < n_individual .or. n_total > m + 2) return

        do si = 1, n_total
            if (si <= n_individual) then
                src = si
            else if (has_sl_agg .and. si == n_individual + 1) then
                src = m + 1
            else
                src = m + 2
            end if

            shifted_rhs = bs_temp(src)
            do sj = 1, qdim
                shifted_rhs = shifted_rhs - gammaxs_temp(src, sj) * theta_offset(sj)
            end do

            if (src == m + 1) then
                if (shifted_rhs > 0.0d0) return
                do sj = 1, qdim
                    gammaxs(si, sj) = -gammaxs_temp(src, sj)
                end do
                bs(si) = -shifted_rhs
                IBs(si) = qdim + n_total + si
            else if (src == m + 2) then
                if (shifted_rhs < 0.0d0) return
                do sj = 1, qdim
                    gammaxs(si, sj) = gammaxs_temp(src, sj)
                end do
                bs(si) = shifted_rhs
                IBs(si) = qdim + si
            else if (shifted_rhs < 0.0d0) then
                do sj = 1, qdim
                    gammaxs(si, sj) = -gammaxs_temp(src, sj)
                end do
                bs(si) = -shifted_rhs
                IBs(si) = qdim + n_total + si
            else
                do sj = 1, qdim
                    gammaxs(si, sj) = gammaxs_temp(src, sj)
                end do
                bs(si) = shifted_rhs
                IBs(si) = qdim + si
            end if
            freevarrow(si) = (si > n_individual)
        end do

        IBs(n_total + 1) = 0
        freevarrow(n_total + 1) = .true.
        bs(n_total + 1) = 0.0d0
        do sj = 1, qdim
            gammaxs(n_total + 1, sj) = 0.0d0
            do si = 1, n_total
                if (IBs(si) > qdim .and. IBs(si) <= qdim + n_total) then
                    gammaxs(n_total + 1, sj) = gammaxs(n_total + 1, sj) - &
                        tau * ws(si) * gammaxs(si, sj)
                else if (IBs(si) > qdim + n_total) then
                    gammaxs(n_total + 1, sj) = gammaxs(n_total + 1, sj) - &
                        (1.0d0 - tau) * ws(si) * gammaxs(si, sj)
                end if
            end do
        end do

        do si = 1, qdim
            r1(si) = si
            r2(si) = 0
        end do
        rr = 0.0d0
        remaining = maxit - iter
        if (remaining <= 0) return
        call run_simplex_full_tvcqr(gammaxs, bs, IBs, freevarrow, r1, r2, rr, ws, &
                                    m + 3, n_total, qdim, tol, remaining, bland, &
                                    iter_fresh, no_pivot_fresh, converged_fresh)
        iter = iter + iter_fresh
        total_simplex_iterations = total_simplex_iterations + iter_fresh
        if (no_pivot_fresh .or. (.not. converged_fresh)) return

        delta_est = 0.0d0
        do si = 1, n_total
            if (IBs(si) >= 1 .and. IBs(si) <= qdim) then
                delta_est(IBs(si)) = bs(si)
            end if
        end do
        do si = 1, qdim
            estimate(si) = theta_offset(si) + delta_est(si)
            reduced_idx = r1(si) - qdim
            if (reduced_idx < 1 .or. reduced_idx > n_individual) return
            H_indices(si) = idx_not_jl_or_jh(reduced_idx)
        end do
        if (.not. tvcqr_H_basis_valid(H_indices, A, m, qdim)) return

        call reorder_fresh_tableau(theta_offset, n_individual, n_total, reorder_ok)
        if (.not. reorder_ok) return
        success = .true.
    end subroutine run_fresh_initialization

    subroutine run_full_active_recovery(success)
        implicit none
        logical, intent(out) :: success
        double precision :: zero_offset(2*(nvar+1))
        integer :: fi, fj

        success = .false.
        zero_offset = 0.0d0
        sl = .false.
        sh = .false.
        has_sl_agg = .false.
        has_sh_agg = .false.
        n_sl = 0
        n_sh = 0
        ms_org = m
        ms = m
        do fi = 1, m
            idx_not_jl_or_jh(fi) = fi
            not_jl_or_jh(fi) = .true.
            ws(fi) = w(fi)
            do fj = 1, 2*(nvar+1)
                gammaxs_temp(fi, fj) = A(fi, fj)
            end do
            bs_temp(fi) = y(fi)
        end do
        call run_fresh_initialization(zero_offset, m, m, success)
    end subroutine run_full_active_recovery

    subroutine run_simplex_full_tvcqr(gx, bv, IBv, fvr, r1v, r2v, rrv, wv, &
                                      ldgx, mv, p, tl, mxit, bld, iters, no_pivot, converged)
        implicit none
        integer, intent(in) :: ldgx, mv, p, mxit
        double precision, intent(inout) :: gx(ldgx, p), bv(mv+1)
        integer, intent(inout) :: IBv(mv+1), r1v(p), r2v(p)
        logical, intent(inout) :: fvr(mv+1)
        double precision, intent(inout) :: rrv(2, p)
        double precision, intent(in) :: wv(mv), tl
        logical, intent(in) :: bld
        integer, intent(out) :: iters
        logical, intent(out) :: no_pivot, converged
        double precision :: yyv(mv+1), eev(mv+1), kval(mv+1)
        double precision :: rrlv, min_kv, pivot_val
        integer :: ii, jj, kk, enter_col, enter_side, enter_var, idx_offset

        iters = 0
        no_pivot = .false.
        converged = .false.
        do while (iters < mxit)
            do ii = 1, p
                rrv(1, ii) = gx(mv+1, ii)
                if (r2v(ii) /= 0) then
                    idx_offset = r1v(ii) - p
                    if (idx_offset >= 1 .and. idx_offset <= mv) then
                        rrv(2, ii) = wv(idx_offset) - rrv(1, ii)
                    else
                        rrv(2, ii) = -rrv(1, ii)
                    end if
                else
                    rrv(2, ii) = 0.0d0
                    rrv(1, ii) = -abs(rrv(1, ii))
                end if
            end do

            rrlv = minval(rrv)
            if (rrlv >= -tl) then
                converged = .true.
                exit
            end if

            enter_col = 0
            enter_side = 0
            enter_var = huge(1)
            if (bld) then
                do jj = 1, p
                    if (rrv(1, jj) < -tl .and. r1v(jj) < enter_var) then
                        enter_var = r1v(jj)
                        enter_col = jj
                        enter_side = 1
                    end if
                end do
                if (enter_col == 0) then
                    do jj = 1, p
                        if (rrv(2, jj) < -tl .and. r2v(jj) < enter_var) then
                            enter_var = r2v(jj)
                            enter_col = jj
                            enter_side = 2
                        end if
                    end do
                end if
            else
                do jj = 1, p
                    do ii = 1, 2
                        if (abs(rrv(ii, jj) - rrlv) < tl) then
                            enter_col = jj
                            enter_side = ii
                            if (ii == 1) then
                                enter_var = r1v(jj)
                            else
                                enter_var = r2v(jj)
                            end if
                            exit
                        end if
                    end do
                    if (enter_col /= 0) exit
                end do
            end if
            if (enter_col == 0) then
                no_pivot = .true.
                exit
            end if

            if (r2v(enter_col) /= 0) then
                if (enter_side == 1) then
                    yyv = gx(1:mv+1, enter_col)
                else
                    yyv = -gx(1:mv+1, enter_col)
                end if
                min_kv = huge(1.0d0)
                kk = 0
                do ii = 1, mv+1
                    if (yyv(ii) > tl .and. .not. fvr(ii)) then
                        kval(ii) = bv(ii) / yyv(ii)
                        if (kval(ii) < min_kv - tl) then
                            min_kv = kval(ii)
                            kk = ii
                        else if (abs(kval(ii) - min_kv) < tl .and. bld) then
                            if (kk == 0 .or. IBv(ii) < IBv(kk)) kk = ii
                        end if
                    end if
                end do
                if (kk == 0) then
                    no_pivot = .true.
                    exit
                end if
                if (enter_side == 2) then
                    idx_offset = r1v(enter_col) - p
                    if (idx_offset >= 1 .and. idx_offset <= mv) then
                        yyv(mv+1) = yyv(mv+1) + wv(idx_offset)
                    end if
                end if
            else
                yyv = gx(1:mv+1, enter_col)
                min_kv = huge(1.0d0)
                kk = 0
                if (yyv(mv+1) < 0.0d0) then
                    do ii = 1, mv+1
                        if (yyv(ii) > tl .and. .not. fvr(ii)) then
                            kval(ii) = bv(ii) / yyv(ii)
                            if (kval(ii) < min_kv - tl) then
                                min_kv = kval(ii)
                                kk = ii
                            else if (abs(kval(ii) - min_kv) < tl .and. bld) then
                                if (kk == 0 .or. IBv(ii) < IBv(kk)) kk = ii
                            end if
                        end if
                    end do
                else
                    do ii = 1, mv+1
                        if (yyv(ii) < -tl .and. .not. fvr(ii)) then
                            kval(ii) = -bv(ii) / yyv(ii)
                            if (kval(ii) < min_kv - tl) then
                                min_kv = kval(ii)
                                kk = ii
                            else if (abs(kval(ii) - min_kv) < tl .and. bld) then
                                if (kk == 0 .or. IBv(ii) < IBv(kk)) kk = ii
                            end if
                        end if
                    end do
                end if
                if (kk == 0) then
                    no_pivot = .true.
                    exit
                end if
                fvr(kk) = .true.
            end if

            do ii = 1, mv+1
                if (ii == kk) then
                    eev(ii) = 1.0d0 - 1.0d0 / yyv(kk)
                else
                    eev(ii) = yyv(ii) / yyv(kk)
                end if
            end do

            if (IBv(kk) <= p + mv) then
                gx(1:mv+1, enter_col) = 0.0d0
                gx(kk, enter_col) = 1.0d0
                r1v(enter_col) = IBv(kk)
                r2v(enter_col) = IBv(kk) + mv
            else
                gx(1:mv+1, enter_col) = 0.0d0
                gx(kk, enter_col) = -1.0d0
                idx_offset = IBv(kk) - p - mv
                if (idx_offset >= 1 .and. idx_offset <= mv) then
                    gx(mv+1, enter_col) = wv(idx_offset)
                end if
                r1v(enter_col) = IBv(kk) - mv
                r2v(enter_col) = IBv(kk)
            end if

            do jj = 1, p
                pivot_val = gx(kk, jj)
                do ii = 1, mv+1
                    gx(ii, jj) = gx(ii, jj) - eev(ii) * pivot_val
                end do
            end do
            pivot_val = bv(kk)
            do ii = 1, mv+1
                bv(ii) = bv(ii) - eev(ii) * pivot_val
            end do
            IBv(kk) = enter_var
            iters = iters + 1
        end do
    end subroutine run_simplex_full_tvcqr

    subroutine handle_current_bad_signs()
        implicit none
        integer :: bi

        if (bad_signs > int(0.1d0 * dble(ms))) then
            mmm_thresh = 2.0d0 * mmm_thresh
            not_new_sl_sh = .true.
        else
            do bi = 1, m
                if ((r(bi) <= 0.0d0) .and. sh(bi)) then
                    sh(bi) = .false.
                end if
                if ((r(bi) >= 0.0d0) .and. sl(bi)) then
                    sl(bi) = .false.
                end if
            end do
            not_new_sl_sh = .false.
        end if
    end subroutine handle_current_bad_signs

    subroutine set_failure(code, eval_idx)
        implicit none
        integer, intent(in) :: code, eval_idx

        ierr = code
        if (eval_idx >= 1 .and. eval_idx <= m) then
            failed_eval = eval_idx
            H_seq(eval_idx, :) = 0
        else
            failed_eval = 1
        end if
    end subroutine set_failure

    logical function certify_tvcqr_candidate(H_idx, r_vec, A_mat, m_loc, p, res_tol_loc)
        implicit none
        integer, intent(in) :: m_loc, p
        integer, intent(in) :: H_idx(p)
        double precision, intent(in) :: r_vec(m_loc), A_mat(m_loc, p), res_tol_loc
        integer :: i, j, ierr_local
        double precision :: AH(p, p), AH_inv(p, p)

        certify_tvcqr_candidate = .true.

        do i = 1, p
            if (H_idx(i) < 1 .or. H_idx(i) > m_loc) then
                certify_tvcqr_candidate = .false.
                return
            end if
            do j = i + 1, p
                if (H_idx(i) == H_idx(j)) then
                    certify_tvcqr_candidate = .false.
                    return
                end if
            end do
            if (abs(r_vec(H_idx(i))) > res_tol_loc) then
                certify_tvcqr_candidate = .false.
                return
            end if
        end do

        do i = 1, p
            do j = 1, p
                AH(i, j) = A_mat(H_idx(i), j)
            end do
        end do
        call matrix_inverse_2p(AH, AH_inv, p, ierr_local)
        if (ierr_local /= 0) then
            certify_tvcqr_candidate = .false.
            return
        end if
    end function certify_tvcqr_candidate

    logical function tvcqr_H_basis_valid(H_idx, A_mat, m_loc, p)
        implicit none
        integer, intent(in) :: m_loc, p
        integer, intent(in) :: H_idx(p)
        double precision, intent(in) :: A_mat(m_loc, p)
        integer :: i, j, ierr_local
        double precision :: AH(p, p), AH_inv(p, p)

        tvcqr_H_basis_valid = .true.

        do i = 1, p
            if (H_idx(i) < 1 .or. H_idx(i) > m_loc) then
                tvcqr_H_basis_valid = .false.
                return
            end if
            do j = i + 1, p
                if (H_idx(i) == H_idx(j)) then
                    tvcqr_H_basis_valid = .false.
                    return
                end if
            end do
        end do

        do i = 1, p
            do j = 1, p
                AH(i, j) = A_mat(H_idx(i), j)
            end do
        end do
        call matrix_inverse_2p(AH, AH_inv, p, ierr_local)
        if (ierr_local /= 0) then
            tvcqr_H_basis_valid = .false.
        end if
    end function tvcqr_H_basis_valid

    subroutine same_h_refit_tvcqr(H_idx, A_mat, y_vec, m_loc, p, res_tol_loc, &
                                  estimate_out, r_out, recovered)
        implicit none
        integer, intent(in) :: m_loc, p
        integer, intent(in) :: H_idx(p)
        double precision, intent(in) :: A_mat(m_loc, p), y_vec(m_loc), res_tol_loc
        double precision, intent(out) :: estimate_out(p), r_out(m_loc)
        logical, intent(out) :: recovered
        integer :: i, j, ierr_local
        double precision :: AH(p, p), AH_inv(p, p), yH(p)

        recovered = .false.
        estimate_out = 0.0d0
        r_out = 0.0d0

        if (.not. tvcqr_H_basis_valid(H_idx, A_mat, m_loc, p)) then
            return
        end if

        do i = 1, p
            yH(i) = y_vec(H_idx(i))
            do j = 1, p
                AH(i, j) = A_mat(H_idx(i), j)
            end do
        end do
        call matrix_inverse_2p(AH, AH_inv, p, ierr_local)
        if (ierr_local /= 0) then
            return
        end if

        do i = 1, p
            estimate_out(i) = 0.0d0
            do j = 1, p
                estimate_out(i) = estimate_out(i) + AH_inv(i, j) * yH(j)
            end do
        end do

        do i = 1, m_loc
            r_out(i) = y_vec(i)
            do j = 1, p
                r_out(i) = r_out(i) - A_mat(i, j) * estimate_out(j)
            end do
        end do

        recovered = certify_tvcqr_candidate(H_idx, r_out, A_mat, m_loc, p, res_tol_loc)
    end subroutine same_h_refit_tvcqr

    ! Helper to find location of minimum in 2D array
    subroutine minloc2d(arr, min_val, row, col)
        double precision, intent(in) :: arr(:,:)
        double precision, intent(in) :: min_val
        integer, intent(out) :: row, col
        integer :: i, j

        do j = 1, size(arr, 2)
            do i = 1, size(arr, 1)
                if (abs(arr(i,j) - min_val) < 1.0d-14) then
                    row = i
                    col = j
                    return
                end if
            end do
        end do
        row = 0
        col = 0
    end subroutine minloc2d
    ! Matrix inversion subroutine using LAPACK
    subroutine matrix_inverse_2p(A_in, A_inv, n, ierr)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: A_in(n, n)
        double precision, intent(out) :: A_inv(n, n)
        integer, intent(out) :: ierr

        ! LAPACK workspace variables
        integer :: lwork, info
        integer :: ipiv(n)
        double precision, allocatable :: work(:)

        ierr = 0

        ! Copy input matrix (LAPACK overwrites)
        A_inv = A_in

        ! LU factorization
        call DGETRF(n, n, A_inv, n, ipiv, info)
        if (info /= 0) then
            ierr = info
            return
        end if

        ! Query optimal workspace
        allocate(work(1))
        lwork = -1
        call DGETRI(n, A_inv, n, ipiv, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! Compute inverse
        call DGETRI(n, A_inv, n, ipiv, work, lwork, info)
        if (info /= 0) then
            ierr = info
            return
        end if

        deallocate(work)

    end subroutine matrix_inverse_2p


    recursive subroutine quicksort_real(arr, left, right)
        implicit none
        double precision, intent(inout) :: arr(:)
        integer, intent(in) :: left, right
        integer :: i, j
        double precision :: pivot, temp

        if (left >= right) return
        if (right - left <= 16) then
            call insertion_sort_real(arr, left, right)
            return
        end if

        pivot = arr((left + right) / 2)
        i = left
        j = right

        do
            do while (arr(i) < pivot)
                i = i + 1
            end do
            do while (arr(j) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                temp = arr(i)
                arr(i) = arr(j)
                arr(j) = temp
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do

        if (left < j) call quicksort_real(arr, left, j)
        if (i < right) call quicksort_real(arr, i, right)
    end subroutine quicksort_real

    subroutine insertion_sort_real(arr, left, right)
        implicit none
        double precision, intent(inout) :: arr(:)
        integer, intent(in) :: left, right
        integer :: i, j
        double precision :: value

        do i = left + 1, right
            value = arr(i)
            j = i - 1
            do while (j >= left)
                if (arr(j) <= value) exit
                arr(j + 1) = arr(j)
                j = j - 1
            end do
            arr(j + 1) = value
        end do
    end subroutine insertion_sort_real

    double precision function median_value(arr, n)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: arr(n)
        double precision :: sorted(n)

        sorted = arr
        call quicksort_real(sorted, 1, n)

        if (mod(n, 2) == 0) then
            median_value = 0.5d0 * (sorted(n/2) + sorted(n/2 + 1))
        else
            median_value = sorted((n+1)/2)
        end if
    end function median_value


end subroutine tvcqr_seq_ppro_fortran
