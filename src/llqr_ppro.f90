! Local Linear Quantile Regression - PPRO Algorithm (M Threshold Warm Start)
! Fortran implementation of llqr_tau_seq_ppro.R
!
! Helper function: 2x2 matrix inverse
subroutine inv22(mat, inv_mat, success)
    implicit none
    double precision, intent(in) :: mat(2, 2)
    double precision, intent(out) :: inv_mat(2, 2)
    logical, intent(out) :: success
    double precision :: det

    det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)

    if (abs(det) < 1.0d-15) then
        success = .false.
        inv_mat = 0.0d0
        return
    end if

    ! Match R's inv22: [d, -b; -c, a] / det
    ! For mat = [[a,b],[c,d]], inv = [[d,-b],[-c,a]] / det
    inv_mat(1,1) = mat(2,2) / det    ! d
    inv_mat(1,2) = -mat(1,2) / det   ! -b
    inv_mat(2,1) = -mat(2,1) / det   ! -c
    inv_mat(2,2) = mat(1,1) / det    ! a
    success = .true.
end subroutine inv22

! Quickselect helper for median computation (average O(n))
subroutine quickselect_inplace(a, n, k, kth_val)
    implicit none
    integer, intent(in) :: n, k
    double precision, intent(inout) :: a(n)
    double precision, intent(out) :: kth_val
    integer :: left, right, i, store, pivot_idx
    double precision :: pivot_val, tmp

    left = 1
    right = n
    do
        if (left == right) then
            kth_val = a(left)
            return
        end if

        pivot_idx = (left + right) / 2
        pivot_val = a(pivot_idx)

        tmp = a(pivot_idx)
        a(pivot_idx) = a(right)
        a(right) = tmp

        store = left
        do i = left, right - 1
            if (a(i) < pivot_val) then
                tmp = a(store)
                a(store) = a(i)
                a(i) = tmp
                store = store + 1
            end if
        end do

        tmp = a(right)
        a(right) = a(store)
        a(store) = tmp

        if (k == store) then
            kth_val = a(store)
            return
        else if (k < store) then
            right = store - 1
        else
            left = store + 1
        end if
    end do
end subroutine quickselect_inplace

! Helper function: median of absolute values (matches R's median(abs(r)))
function median_abs(arr, n) result(median_val)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: arr(n)
    double precision :: median_val
    double precision :: work(n), v1, v2
    integer :: i, k1, k2

    do i = 1, n
        work(i) = abs(arr(i))
    end do

    k1 = (n + 1) / 2
    k2 = (n + 2) / 2
    call quickselect_inplace(work, n, k1, v1)
    if (k1 == k2) then
        median_val = v1
    else
        do i = 1, n
            work(i) = abs(arr(i))
        end do
        call quickselect_inplace(work, n, k2, v2)
        median_val = 0.5d0 * (v1 + v2)
    end if
end function median_abs

! Helper function: max of array
function max_array(arr, n) result(maxval)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: arr(n)
    double precision :: maxval
    integer :: i

    maxval = arr(1)
    do i = 2, n
        if (arr(i) > maxval) maxval = arr(i)
    end do
end function max_array

! Main PPRO subroutine
subroutine llqr_ppro_fortran(x, y, z, m, nvar, rounds, tau, h, tol, maxit, &
                             Mm_factor, case_int, bland_int, min_subsample_size_in, &
                             ll_est, d_ll_est, H_mat, full_active_recovery_count, &
                             ierr, failed_eval, always_same_h_refit_int, &
                             threshold_lower_bound_int, threshold_scale_mode_int)

    implicit none

    ! Input arguments
    integer, intent(in) :: m, nvar, rounds, maxit, case_int, bland_int
    integer, intent(in) :: min_subsample_size_in, always_same_h_refit_int
    integer, intent(in) :: threshold_lower_bound_int, threshold_scale_mode_int
    double precision, intent(in) :: x(m), y(m), z(rounds), tau, tol, Mm_factor
    double precision, intent(inout) :: h

    ! Output arguments
    double precision, intent(out) :: ll_est(rounds)
    double precision, intent(out) :: d_ll_est(rounds)
    integer, intent(out) :: H_mat(rounds, nvar+1)
    integer, intent(out) :: full_active_recovery_count
    integer, intent(out) :: ierr
    integer, intent(out) :: failed_eval

    ! Local variables for full problem (round 1)
    double precision :: A(m, nvar+1)     ! Design matrix [1, x]
    double precision :: w(m)              ! Kernel weights
    double precision :: eva_z(m)          ! z - x for kernel
    double precision :: gammax(m+1, nvar+1)
    double precision :: b(m+1)
    integer :: IB(m+1)
    logical :: freevarrow(m+1)
    integer :: r1(nvar+1), r2(nvar+1)
    double precision :: rr(2, nvar+1)

    ! Working variables
    double precision :: yy(m+1), ee(m+1), k_vals(m+1)
    double precision :: u(m), v(m), estimate(nvar+1)
    integer :: i, j, k, rd, iter, t_rr, tsep, t, jj
    integer :: iter_total, iter_attempt, remaining
    integer :: idx_r1_minus_offset
    double precision :: rrl, min_k, pi, pivot_row_value
    logical :: bland, always_same_h_refit, threshold_lower_bound

    ! PPRO-specific variables
    double precision :: mm, mmm, M_threshold, residual_scale, threshold_scale
    logical :: sl(m), sh(m), not_jl_or_jh(m), active(m)
    logical :: has_sl_agg, has_sh_agg
    integer :: ms  ! subsample size
    integer :: H_prev(nvar+1)
    integer :: n_bad_signs
    logical :: not_optimal, not_new_sl_sh  ! Bad signs loop control
    logical :: force_full_sample, no_pivot_flag, accept_subsample
    logical :: simplex_converged
    integer :: empty_pivot_count, max_empty_pivot_retries
    double precision :: res_tol

    ! Subsample arrays (will be allocated dynamically in concept, but use max size)
    double precision :: gammaxs(m+2, nvar+1)  ! max possible size
    double precision :: bs(m+2)
    double precision :: ws(m+1)
    integer :: IBs(m+2)
    logical :: freevarrows(m+3)
    integer :: r1s(nvar+1), r2s(nvar+1)

    ! Temporary storage
    double precision :: gammaxs_temp(m+2, nvar+1)
    double precision :: bs_temp(m+2)
    integer :: idx_not_jl_or_jh(m)
    integer :: n_subsample, n_potential_S
    integer :: min_subsample_size, min_subsample_size_effective
    integer :: n_active, target_min
    double precision :: r(m), r_raw(m)
    double precision :: sum_w_sl, sum_w_sh
    double precision :: ll_candidate, d_ll_candidate
    integer :: H_candidate(nvar+1)
    integer :: h_failure_code
    logical :: h_map_ok
    logical :: cert_reject_return
    logical :: refit_recertified
    logical :: use_independent_init, H_prev_valid
    logical :: shifted_init_success, full_recovery_success
    integer :: independent_trigger
    double precision :: theta_prev(nvar+1)

    ! Z-sorting variables (CRITICAL FIX: match R's z-sorting behavior)
    double precision :: z_sorted(rounds)
    integer :: z_order(rounds)
    double precision :: ll_est_sorted(rounds), d_ll_est_sorted(rounds)
    double precision :: residual_prev(m)
    integer :: H_mat_sorted(rounds, nvar+1)

    ! Warm start variables (for rd==2)
    double precision :: xh(nvar+1, nvar+1)
    double precision :: xhinv(nvar+1, nvar+1)
    logical :: inv_success
    integer :: H_subsample(nvar+1)
    integer :: idpos(m), idneg(m)
    integer :: n_idpos, n_idneg, n_hbar, n_hbar_core
    integer :: Hbar(m)
    double precision :: P(m)
    integer :: u_in_IBs(m), v_in_IBs(m)
    double precision :: Pxhbar(m, nvar+1)
    double precision :: Pxhbarxhinv(m, nvar+1)
    double precision :: bs_hbar(m)
    double precision :: lambda(m)
    double precision :: obj_row(nvar+1)
    integer :: r_idx
    logical :: is_in_H
    double precision :: u_subsample(m), v_subsample(m)
    ! Temporary arrays for simplex with correct dimensions
    double precision :: gammaxs_simplex(m+1, nvar+1)
    double precision :: bs_simplex(m+1)

    ! rd>2 incremental update variables (matching R's approach)
    double precision :: xhinv_stored(nvar+1, nvar+1)
    double precision :: bs_stored(nvar+1)
    double precision :: gammaxs_pos(m, nvar+1)  ! Stored positive Hbar rows
    double precision :: bs_pos(m)
    integer :: idx_Hbar_pos(m), n_pos_prev
    double precision :: gammaxs_neg(m, nvar+1)  ! Stored negative Hbar rows
    double precision :: bs_neg(m)
    integer :: idx_Hbar_neg(m), n_neg_prev
    integer :: idx_Hbar_pos2(m), idx_Hbar_neg2(m)
    logical :: matched_rows_pos(m), matched_rows_neg(m)
    integer :: n_matched_pos, n_matched_neg, n_unmatched_pos, n_unmatched_neg
    integer :: matched_idx_pos(m), matched_idx_neg(m)
    integer :: unmatched_idx_pos(m), unmatched_idx_neg(m)
    integer :: row_map_pos(m), row_map_neg(m)
    double precision :: Pxhbarxhinv_new(m, nvar+1)
    double precision :: temp_vec(nvar+1)
    integer :: ms_org, curr_idx, prev_idx, match_pos
    integer :: n_xhinv_hbar, ms_current, n_lambda  ! For rd>2 tableau assembly

    ! Helper function declarations
    double precision :: median_abs, max_array

    ! Constants
    pi = 4.0d0 * atan(1.0d0)
    bland = (bland_int /= 0)
    always_same_h_refit = (always_same_h_refit_int /= 0)
    threshold_lower_bound = (threshold_lower_bound_int /= 0)
    res_tol = 1.0d-6
    max_empty_pivot_retries = 3
    ierr = 0
    failed_eval = 0
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

    if (min_subsample_size_in > 0) then
        min_subsample_size_effective = min_subsample_size_in
    else
        min_subsample_size_effective = max(5*(nvar + 1), ceiling(0.2d0 * dble(m)))
    end if

    ! Set default bandwidth
    if (h <= 0.0d0) then
        h = dble(m)**(-0.2d0)
    end if

    ! ============================================================
    ! CRITICAL FIX: Sort z array (matching R behavior!)
    ! R code: original_order <- order(z); z <- z[original_order]
    ! This is ESSENTIAL for PPRO's warm-start to work correctly!
    ! ============================================================
    do i = 1, rounds
        z_sorted(i) = z(i)
        z_order(i) = i
    end do

    ! Simple bubble sort (sufficient for typical evaluation points)
    do i = 1, rounds-1
        do j = i+1, rounds
            if (z_sorted(i) > z_sorted(j)) then
                ! Swap z values
                min_k = z_sorted(i)
                z_sorted(i) = z_sorted(j)
                z_sorted(j) = min_k
                ! Swap order indices
                k = z_order(i)
                z_order(i) = z_order(j)
                z_order(j) = k
            end if
        end do
    end do

    ! ============================================================
    ! Build design matrix A = [1, x] (constant across eval points)
    ! ============================================================
    do i = 1, m
        A(i, 1) = 1.0d0
        A(i, 2) = x(i)
    end do

    ! Match R llqr_seq_ppro threshold order by data-generation case
    if (case_int == 1) then
        mm = log(log(dble(m))) / sqrt(log(dble(m)))
    else if (case_int == 2) then
        ! Use the actual bandwidth.  This reduces to
        ! sqrt(log(m)) * m**(-0.4) when h = m**(-0.2).
        mm = sqrt(log(dble(m)) / (dble(m) * h))
    else
        ierr = 5
        failed_eval = 1
        return
    end if

    ! ============================================================
    ! ROUND 1: Standard simplex (identical to seq)
    ! ============================================================
    rd = 1

    ! Compute kernel weights for first z (USING SORTED Z!)
    do i = 1, m
        eva_z(i) = z_sorted(rd) - x(i)
        w(i) = llqr_kernel_weight(eva_z(i) / h, case_int, pi)
    end do

    ! Initialize gammax from design matrix
    do i = 1, m
        do j = 1, nvar+1
            gammax(i, j) = A(i, j)
        end do
    end do

    ! Flip signs for negative y
    do i = 1, m
        if (y(i) < 0.0d0) then
            gammax(i, :) = -gammax(i, :)
        end if
    end do

    ! Initialize b vector
    do i = 1, m
        b(i) = abs(y(i))
    end do
    b(m+1) = 0.0d0

    ! Initialize basis IB
    do i = 1, m
        if (y(i) >= 0.0d0) then
            IB(i) = i + nvar + 1  ! u_i
        else
            IB(i) = i + nvar + 1 + m  ! v_i
        end if
    end do
    ! Objective row (m+1) is protected by freevarrow, so IB(m+1) doesn't matter
    ! Set to 0 to avoid conflicts with actual variable indices
    IB(m+1) = 0

    ! Initialize freevarrow
    do i = 1, m
        freevarrow(i) = .false.
    end do
    freevarrow(m+1) = .true.

    ! Initialize r1, r2 (non-basic variables)
    do i = 1, nvar+1
        r1(i) = i
        r2(i) = 0
    end do

    rr = 0.0d0

    ! Compute objective row for round 1
    do j = 1, nvar+1
        gammax(m+1, j) = 0.0d0
        do i = 1, m
            if (IB(i) > nvar + 1 .and. IB(i) <= nvar + 1 + m) then
                ! u_i is basic: use tau * w[i]
                idx_r1_minus_offset = IB(i) - nvar - 1
                if (idx_r1_minus_offset >= 1 .and. idx_r1_minus_offset <= m) then
                    gammax(m+1, j) = gammax(m+1, j) - tau * w(idx_r1_minus_offset) * gammax(i, j)
                end if
            else if (IB(i) > nvar + 1 + m) then
                ! v_i is basic: use (1-tau) * w[i]
                idx_r1_minus_offset = IB(i) - nvar - 1 - m
                if (idx_r1_minus_offset >= 1 .and. idx_r1_minus_offset <= m) then
                    gammax(m+1, j) = gammax(m+1, j) - (1.0d0 - tau) * w(idx_r1_minus_offset) * gammax(i, j)
                end if
            end if
        end do
    end do

    ! Run the full cold-start simplex for round 1.
    call run_simplex_full_llqr(gammax, b, IB, freevarrow, r1, r2, rr, w, m+1, m, nvar, &
                               tau, tol, maxit, bland, iter, no_pivot_flag, simplex_converged)
    if (no_pivot_flag .or. (.not. simplex_converged)) then
        call set_failure(2)
        return
    end if

    ! Extract solution for round 1
    call extract_solution(gammax, b, IB, m, nvar, estimate, u, v, r1)

    ! Store results in SORTED arrays (will unsort at the end)
    ll_est_sorted(rd) = estimate(1) + estimate(2) * z_sorted(rd)
    d_ll_est_sorted(rd) = estimate(2)
    residual_prev = u - v
    H_mat_sorted(rd, :) = r1 - 1 - nvar

    ! BUG FIX: Do NOT sort H_mat for rd=1!
    ! R code keeps H in simplex order: H <- r1 - 1 - nvar (no sort)
    ! The order affects xhinv = solve(A[H,]) and ws[H] in subsequent rounds
    ! Sorting was causing rd=2 to use different H indices than R


    ! ============================================================
    ! ROUNDS 2+: PPRO with M threshold
    ! ============================================================
    do rd = 2, rounds
        ! Initialize bad signs loop control
        not_optimal = .true.
        not_new_sl_sh = .true.
        force_full_sample = .false.
        empty_pivot_count = 0
        mmm = mm
        iter_total = 0
        use_independent_init = .false.
        independent_trigger = 0

        ! Compute kernel weights for current z (USING SORTED Z!)
        n_active = 0
        do i = 1, m
            eva_z(i) = z_sorted(rd) - x(i)
            w(i) = llqr_kernel_weight(eva_z(i) / h, case_int, pi)
            active(i) = (w(i) > 0.0d0)
            if (active(i)) n_active = n_active + 1
        end do

        ! BAD SIGNS OUTER LOOP: Keep trying until solution is good
        attempt_loop: do while (not_optimal)

            ! ========================================================
            ! Step 1: Compute M threshold from previous residuals
            ! ========================================================
            r = residual_prev
            H_prev = H_mat_sorted(rd-1, :)
            if (.not. use_independent_init) then
                call validate_previous_H(H_prev, H_prev_valid)
                if (.not. H_prev_valid) then
                    use_independent_init = .true.
                    independent_trigger = 1
                end if
            end if
            if (not_new_sl_sh) then
                residual_scale = median_abs(r, m)
                if (threshold_lower_bound) then
                    M_threshold = max(Mm_factor * mmm * threshold_scale, 0.1d0 * residual_scale)
                else
                    M_threshold = Mm_factor * mmm * threshold_scale
                end if

                min_subsample_size = min_subsample_size_effective

                target_min = min(min_subsample_size, n_active)
                min_k = 0.0d0
                n_potential_S = 0
                do i = 1, m
                    if (active(i)) then
                        if (abs(r(i)) > min_k) min_k = abs(r(i))
                        if (abs(r(i)) <= M_threshold) n_potential_S = n_potential_S + 1
                    end if
                end do
                do while (n_potential_S < target_min .and. M_threshold < min_k)
                    M_threshold = M_threshold * 1.5d0
                    n_potential_S = 0
                    do i = 1, m
                        if (active(i) .and. abs(r(i)) <= M_threshold) then
                            n_potential_S = n_potential_S + 1
                        end if
                    end do
                end do

                ! Classify observations
                do i = 1, m
                    sl(i) = active(i) .and. r(i) < -M_threshold
                    sh(i) = active(i) .and. r(i) > M_threshold
                end do
            end if

            ! Force H observations into the retained subsample every attempt,
            ! including after few-bad-sign repairs.
            if (force_full_sample) then
                sl = .false.
                sh = .false.
                not_jl_or_jh = active
            else
                do i = 1, m
                    not_jl_or_jh(i) = active(i) .and. (.not. (sl(i) .or. sh(i)))
                end do
            end if

            ! Previous-H rows are zero-cost basis padding when their current
            ! kernel weight is zero; keep them out of active screening counts.
            if (.not. use_independent_init) then
                do i = 1, nvar+1
                    if (sl(H_prev(i)) .or. sh(H_prev(i))) then
                        sl(H_prev(i)) = .false.
                        sh(H_prev(i)) = .false.
                    end if
                    not_jl_or_jh(H_prev(i)) = .true.
                end do
            end if

            sum_w_sl = 0.0d0
            sum_w_sh = 0.0d0
            do i = 1, m
                if (sl(i)) sum_w_sl = sum_w_sl + w(i)
                if (sh(i)) sum_w_sh = sum_w_sh + w(i)
            end do
            has_sl_agg = (sum_w_sl > 0.0d0)
            has_sh_agg = (sum_w_sh > 0.0d0)

        ! Count subsample size
        ms = count(not_jl_or_jh)
        n_subsample = ms

        ! Get indices of subsample observations
        j = 0
        do i = 1, m
            if (not_jl_or_jh(i)) then
                j = j + 1
                idx_not_jl_or_jh(j) = i
            end if
        end do


        ! ========================================================
        ! Step 2: Build subsample problem
        ! ========================================================
        ! NOTE: Subsample observation rows are built ONLY for rd==2 (R line 567)
        ! For rd>2, only aggregates are built; subsample rows are NOT rebuilt!
        ms = n_subsample

        if (rd == 2 .or. use_independent_init) then
            ! Shifted independent initialization also needs raw retained rows.
            do i = 1, ms
                do j = 1, nvar+1
                    gammaxs_temp(i, j) = A(idx_not_jl_or_jh(i), j)
                end do
                bs_temp(i) = y(idx_not_jl_or_jh(i))
            end do
        end if

        ! Copy weights for all rd>=2
        do i = 1, n_subsample
            ws(i) = w(idx_not_jl_or_jh(i))
        end do

        ! Build aggregates for both rd==2 and rd>2 (R lines 572-591)
        if (rd >= 2) then
            ! Aggregate sl observations if any
            ! NOTE: sl is a logical array of size m, so loop over ALL m observations
            ! Store at index m+1 to match R's storage at n+1
            ! R: glob.wx <- colSums(gammaxsl * wsl) where wsl = w[sl]
            ! R: glob.wy <- sum(y[sl] * wsl)
            if (has_sl_agg) then
                do j = 1, nvar+1
                    gammaxs_temp(m+1, j) = 0.0d0
                    do i = 1, m
                        if (sl(i)) then
                            gammaxs_temp(m+1, j) = gammaxs_temp(m+1, j) + A(i, j) * w(i)
                        end if
                    end do
                end do
                bs_temp(m+1) = 0.0d0
                do i = 1, m
                    if (sl(i)) then
                        bs_temp(m+1) = bs_temp(m+1) + y(i) * w(i)
                    end if
                end do
                ms = ms + 1
                ws(ms) = 1.0d0  ! Weight for sl aggregate is 1
            end if

            ! Aggregate sh observations if any
            ! NOTE: sh is a logical array of size m, so loop over ALL m observations
            ! Store at index m+2 to match R's storage at n+2
            ! R: glob.wx <- colSums(gammaxsh * wsh) where wsh = w[sh]
            ! R: glob.wy <- sum(y[sh] * wsh)
            if (has_sh_agg) then
                do j = 1, nvar+1
                    gammaxs_temp(m+2, j) = 0.0d0
                    do i = 1, m
                        if (sh(i)) then
                            gammaxs_temp(m+2, j) = gammaxs_temp(m+2, j) + A(i, j) * w(i)
                        end if
                    end do
                end do
                bs_temp(m+2) = 0.0d0
                do i = 1, m
                    if (sh(i)) then
                        bs_temp(m+2) = bs_temp(m+2) + y(i) * w(i)
                    end if
                end do
                ms = ms + 1
                ws(ms) = 1.0d0  ! Weight for sh aggregate is 1
            end if

        end if

        ! ========================================================
        ! Step 3: Independent or transported-basis initialization
        ! ========================================================
        if (use_independent_init) then
            theta_prev(2) = d_ll_est_sorted(rd - 1)
            theta_prev(1) = ll_est_sorted(rd - 1) - theta_prev(2) * z_sorted(rd - 1)
            call run_shifted_reduced_initialization(shifted_init_success)

            if (.not. shifted_init_success) then
                full_active_recovery_count = full_active_recovery_count + 1
                call run_full_active_recovery(full_recovery_success)
                if (.not. full_recovery_success) then
                    call set_failure(2)
                    return
                end if

                ll_candidate = estimate(1) + estimate(2) * z_sorted(rd)
                d_ll_candidate = estimate(2)
                do i = 1, m
                    r_raw(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
                end do
                if (always_same_h_refit) then
                    call try_same_h_refit(refit_recertified)
                    if (.not. refit_recertified) then
                        do i = 1, m
                            r_raw(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
                        end do
                    end if
                end if
                accept_subsample = certify_llqr_candidate(H_candidate, r_raw, A, m, nvar, res_tol)
                if (.not. accept_subsample) then
                    call set_failure(1)
                    return
                end if

                ll_est_sorted(rd) = ll_candidate
                d_ll_est_sorted(rd) = d_ll_candidate
                H_mat_sorted(rd, :) = H_candidate
                r = r_raw
                do i = 1, nvar + 1
                    r(H_candidate(i)) = 0.0d0
                end do
                residual_prev = r
                call store_current_cache()
                not_optimal = .false.
                cycle attempt_loop
            end if

            ll_candidate = estimate(1) + estimate(2) * z_sorted(rd)
            d_ll_candidate = estimate(2)
            do i = 1, m
                r_raw(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
            end do

            if (always_same_h_refit) then
                call try_same_h_refit(refit_recertified)
                if (.not. refit_recertified) then
                    do i = 1, m
                        r_raw(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
                    end do
                end if
            end if

            n_bad_signs = 0
            do i = 1, m
                if ((sh(i) .and. r_raw(i) <= 0.0d0) .or. (sl(i) .and. r_raw(i) >= 0.0d0)) then
                    n_bad_signs = n_bad_signs + 1
                end if
            end do

            if (n_bad_signs > 0) then
                if (dble(n_bad_signs) > 0.1d0 * dble(ms)) then
                    mmm = mmm * 2.0d0
                    not_new_sl_sh = .true.
                else
                    do i = 1, m
                        if (sh(i) .and. r_raw(i) <= 0.0d0) sh(i) = .false.
                        if (sl(i) .and. r_raw(i) >= 0.0d0) sl(i) = .false.
                    end do
                    not_new_sl_sh = .false.
                end if
                cycle attempt_loop
            end if

            accept_subsample = certify_llqr_candidate(H_candidate, r_raw, A, m, nvar, res_tol)
            if ((.not. accept_subsample) .and. (.not. always_same_h_refit)) then
                call try_same_h_refit(refit_recertified)
                if (refit_recertified) accept_subsample = .true.
            end if
            if (.not. accept_subsample) then
                call handle_certification_reject(cert_reject_return)
                if (cert_reject_return) return
                cycle attempt_loop
            end if

            ll_est_sorted(rd) = ll_candidate
            d_ll_est_sorted(rd) = d_ll_candidate
            H_mat_sorted(rd, :) = H_candidate
            r = r_raw
            do i = 1, nvar + 1
                r(H_candidate(i)) = 0.0d0
            end do
            residual_prev = r
            call store_current_cache()
            not_optimal = .false.
            cycle attempt_loop

        else if (rd == 2) then
            ! Map H from previous round to subsample coordinates
            do i = 1, nvar+1
                H_subsample(i) = 0
                do j = 1, n_subsample
                    if (idx_not_jl_or_jh(j) == H_prev(i)) then
                        H_subsample(i) = j
                        exit
                    end if
                end do
                if (H_subsample(i) == 0) then
                    use_independent_init = .true.
                    independent_trigger = 2
                    not_new_sl_sh = .true.
                    cycle attempt_loop
                end if
            end do

            ! Extract xh = gammaxs[H,]
            do i = 1, nvar+1
                do j = 1, nvar+1
                    xh(i, j) = gammaxs_temp(H_subsample(i), j)
                end do
            end do

            ! Compute inverse
            call inv22(xh, xhinv, inv_success)

            if (.not. inv_success) then
                use_independent_init = .true.
                independent_trigger = 3
                not_new_sl_sh = .true.
                cycle attempt_loop
            end if

            ! Identify positive and negative residuals in subsample (excluding H)
            n_idpos = 0
            n_idneg = 0
            do i = 1, n_subsample
                ! Check if this index is in H
                is_in_H = .false.
                do j = 1, nvar+1
                    if (H_subsample(j) == i) then
                        is_in_H = .true.
                        exit
                    end if
                end do

                if (.not. is_in_H) then
                    ! Get residual from previous round for this subsample observation
                    ! Use the previous accepted residual vector.
                    r_idx = idx_not_jl_or_jh(i)
                    if (residual_prev(r_idx) > 0.0d0) then
                        n_idpos = n_idpos + 1
                        idpos(n_idpos) = i
                    else if (residual_prev(r_idx) < 0.0d0) then
                        n_idneg = n_idneg + 1
                        idneg(n_idneg) = i
                    else
                    end if
                end if
            end do



            ! Build Hbar (concatenate idpos, idneg, and aggregate indices)
            ! R code: Hbar <- c(idpos,idneg,ms-1,ms) when both sl and sh exist
            n_hbar = n_idpos + n_idneg
            do i = 1, n_idpos
                Hbar(i) = idpos(i)
            end do
            do i = 1, n_idneg
                Hbar(n_idpos + i) = idneg(i)
            end do

            ! Add aggregate row indices if they exist
            ! NOTE: In Fortran, Hbar contains indices into gammaxs_temp for accessing original rows
            ! Aggregates are stored at m+1 (sl) and m+2 (sh) in gammaxs_temp
            ! The Fortran code accesses gammaxs_temp(Hbar(i), :) and stores result in gammaxs(nvar+1+i, :)
            ! So Hbar should contain temp array indices, not final tableau indices
            if (has_sl_agg .and. has_sh_agg) then
                ! Both sl and sh exist
                Hbar(n_hbar + 1) = m + 1  ! sl aggregate in temp
                Hbar(n_hbar + 2) = m + 2  ! sh aggregate in temp
                n_hbar = n_hbar + 2
            else if (has_sl_agg) then
                ! Only sl exists
                Hbar(n_hbar + 1) = m + 1  ! sl aggregate in temp
                n_hbar = n_hbar + 1
            else if (has_sh_agg) then
                ! Only sh exists
                Hbar(n_hbar + 1) = m + 2  ! sh aggregate in temp
                n_hbar = n_hbar + 1
            end if

            ! Build P vector (1 for idpos, -1 for idneg, sign for aggregates)
            do i = 1, n_idpos
                P(i) = 1.0d0
            end do
            do i = 1, n_idneg
                P(n_idpos + i) = -1.0d0
            end do

            ! Add P values for aggregates
            if (has_sl_agg .and. has_sh_agg) then
                P(n_idpos + n_idneg + 1) = -1.0d0  ! sl aggregate gets -1
                P(n_idpos + n_idneg + 2) = 1.0d0   ! sh aggregate gets +1
            else if (has_sl_agg) then
                P(n_idpos + n_idneg + 1) = -1.0d0  ! sl aggregate gets -1
            else if (has_sh_agg) then
                P(n_idpos + n_idneg + 1) = 1.0d0   ! sh aggregate gets +1
            end if



            ! Build u_in_IBs and v_in_IBs
            ! R code: u.in.IBs <- idpos + (nvar + 1)
            !         v.in.IBs <- idneg + (nvar + 1) + ms
            ! FIXED: Use ms (with aggregates), not n_subsample!
            do i = 1, n_idpos
                u_in_IBs(i) = idpos(i) + nvar + 1
            end do
            do i = 1, n_idneg
                v_in_IBs(i) = idneg(i) + nvar + 1 + ms
            end do


            ! CRITICAL: Initialize IBs array to avoid garbage values
            ! This prevents segfault from uninitialized memory access
            ! Note: freevarrows is initialized separately below (lines 706-714)
            do i = 1, m+2
                IBs(i) = nvar + 1  ! Safe sentinel value
            end do

            ! Build IBs (initial basis)
            ! R code: IBs <- c(1:(nvar + 1), u.in.IBs, v.in.IBs, nvar+1+2*ms-1, nvar+1+ms)
            ! First nvar+1 entries: 1:(nvar+1) - beta coefficients are basic
            do i = 1, nvar+1
                IBs(i) = i
            end do
            ! Next: u_in_IBs
            do i = 1, n_idpos
                IBs(nvar + 1 + i) = u_in_IBs(i)
            end do
            ! Next: v_in_IBs
            do i = 1, n_idneg
                IBs(nvar + 1 + n_idpos + i) = v_in_IBs(i)
            end do
            ! Next: Add aggregate variable indices if they exist
            ! R code: if both sl and sh exist, add nvar+1+2*ms-1 (v_L for sl) and nvar+1+ms (u_H for sh)
            if (has_sl_agg .and. has_sh_agg) then
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + 2*ms - 1  ! v_L for sl aggregate
                IBs(nvar + 1 + n_idpos + n_idneg + 2) = nvar + 1 + ms        ! u_H for sh aggregate
            else if (has_sl_agg) then
                ! Only sl: R code uses nvar+1+2*ms for v_L
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + 2*ms
            else if (has_sh_agg) then
                ! Only sh: R code uses nvar+1+ms for u_H
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + ms
            end if
            ! REMOVED: Don't add dummy for objective row - R doesn't do this
            ! IBs should have exactly ms elements, matching R


            ! Initialize freevarrows to match the R construction:
            ! c(rep(TRUE, nvar+1), rep(FALSE, length(u.in.IBs)), rep(FALSE, length(v.in.IBs)),
            !   [TRUE for each aggregate row], TRUE for the objective row)
            n_hbar_core = n_idpos + n_idneg
            do i = 1, ms + 1
                freevarrows(i) = .false.
            end do
            do i = 1, nvar+1
                freevarrows(i) = .true.
            end do
            do i = nvar + 2, nvar + 1 + n_hbar_core
                freevarrows(i) = .false.
            end do
            if (has_sl_agg .and. has_sh_agg) then
                freevarrows(nvar + 1 + n_hbar_core + 1) = .true.  ! sl aggregate row
                freevarrows(nvar + 1 + n_hbar_core + 2) = .true.  ! sh aggregate row
            else if (has_sl_agg .or. has_sh_agg) then
                freevarrows(nvar + 1 + n_hbar_core + 1) = .true.  ! single aggregate row
            end if
            freevarrows(ms + 1) = .true.  ! objective row

            ! Initialize r1s and r2s for the reduced problem
            ! Non-basic variables are H observations
            ! NOTE: r2 uses ms (subsample size WITH aggregates), not n_subsample
            do i = 1, nvar+1
                r1s(i) = H_subsample(i) + nvar + 1
                r2s(i) = r1s(i) + ms
            end do

            ! Build warm start tableau
            ! First, copy non-H subsample rows from gammaxs_temp to gammaxs
            ! R code: gammaxs <- gammaxs.temp[c(1:(ms-2), n+1, n+2),]
            ! Rows 1:(ms-2) are non-H subsample observations
            ! NOTE: In Fortran, subsample rows are in positions 1:n_subsample of gammaxs_temp,
            ! but aggregates are at m+1 and m+2, not at n_subsample+1 and n_subsample+2!
            do i = 1, n_subsample
                do j = 1, nvar+1
                    gammaxs(i, j) = gammaxs_temp(i, j)
                end do
            end do

            ! Then append aggregate rows from m+1 and m+2 if they exist
            if (has_sl_agg) then
                do j = 1, nvar+1
                    gammaxs(n_subsample + 1, j) = gammaxs_temp(m+1, j)
                end do
            end if
            if (has_sh_agg) then
                k = n_subsample + 1
                if (has_sl_agg) k = k + 1  ! If sl exists, sh goes to next position
                do j = 1, nvar+1
                    gammaxs(k, j) = gammaxs_temp(m+2, j)
                end do
            end if

            ! Then overwrite first nvar+1 rows with xhinv
            do i = 1, nvar+1
                do j = 1, nvar+1
                    gammaxs(i, j) = xhinv(i, j)
                end do
            end do

            ! Compute bs[1:(nvar+1)] = xhinv %*% bs_temp[H]
            do i = 1, nvar+1
                bs(i) = 0.0d0
                do j = 1, nvar+1
                    bs(i) = bs(i) + xhinv(i, j) * bs_temp(H_subsample(j))
                end do
            end do



            ! Hbar rows: -P * gammaxs[Hbar,] * xhinv
            if (n_hbar > 0) then
                do i = 1, nvar+1
                end do

                ! Compute Pxhbarxhinv and bs row by row, storing each immediately
                ! IMPORTANT: R computes Pxhbarxhinv WITHOUT leading negative sign!
                ! R code:
                !   Pxhbar <- gammaxs[Hbar, ] * P  (element-wise, NO negative!)
                !   Pxhbarxhinv <- Pxhbar %*% xhinv
                !   gammaxs <- rbind(xhinv, - Pxhbarxhinv)  (negative added when storing to tableau)
                ! For objective row: obj_row <- tau * ws[H] + t(lambda) %*% Pxhbarxhinv (uses Pxhbarxhinv WITHOUT negative)
                do i = 1, n_hbar
                    ! Step 1: Compute P[i] * gammaxs_ORIGINAL[Hbar[i],]
                    ! FIXED: Removed incorrect leading negative sign
                    ! IMPORTANT: Must use gammaxs_temp (original design matrix), NOT gammaxs!
                    ! Because gammaxs[1:2,] has been overwritten with xhinv
                    do j = 1, nvar+1
                        Pxhbar(i, j) = P(i) * gammaxs_temp(Hbar(i), j)
                    end do

                    ! Step 2: Multiply by xhinv to get Pxhbarxhinv[i,]
                    do j = 1, nvar+1
                        Pxhbarxhinv(i, j) = 0.0d0
                        do k = 1, nvar+1
                            Pxhbarxhinv(i, j) = Pxhbarxhinv(i, j) + Pxhbar(i, k) * xhinv(k, j)
                        end do
                    end do

                    ! Step 3: IMMEDIATELY store into gammaxs[nvar+1+i,]
                    ! R stores -Pxhbarxhinv in the tableau: gammaxs <- rbind(xhinv, - Pxhbarxhinv)
                    do j = 1, nvar+1
                        gammaxs(nvar + 1 + i, j) = -Pxhbarxhinv(i, j)
                    end do

                    ! Step 4: Compute bs[nvar+1+i]
                    ! R formula: bs <- c(xhinv %*% bs[H], - Pxhbarxhinv %*% bs[H] + bs[Hbar] * P)
                    ! FIXED: Use Pxhbarxhinv, not Pxhbar!
                    ! bs[nvar+1+i] = -Pxhbarxhinv[i,] @ bs_temp[H] + bs_temp[Hbar[i]] * P[i]
                    bs_hbar(i) = 0.0d0
                    do j = 1, nvar+1
                        bs_hbar(i) = bs_hbar(i) - Pxhbarxhinv(i, j) * bs_temp(H_subsample(j))
                    end do
                    bs_hbar(i) = bs_hbar(i) + bs_temp(Hbar(i)) * P(i)

                    ! Step 5: IMMEDIATELY store bs[nvar+1+i]
                    bs(nvar + 1 + i) = bs_hbar(i)
                end do

                ! Compute objective row
                ! lambda = [tau * ws[idpos], (1-tau) * ws[idneg], (1-tau), tau]
                ! where the last two elements are for sl and sh aggregates
                do i = 1, n_idpos
                    lambda(i) = tau * ws(idpos(i))
                end do
                do i = 1, n_idneg
                    lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                end do
                ! Add aggregate elements to lambda
                ! NOTE: Hbar is [idpos, idneg, m+1, m+2] so aggregates are at positions n_idpos+n_idneg+1 and +2
                if (has_sl_agg) then
                    lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau  ! sl aggregate: 1-tau
                end if
                if (has_sh_agg) then
                    k = n_idpos + n_idneg + 1
                    if (has_sl_agg) k = k + 1
                    lambda(k) = tau  ! sh aggregate: tau
                end if

                ! obj_row = tau * ws[H] + t(lambda) %*% Pxhbarxhinv
                ! R code: obj_row <- tau * ws[H] + crossprod(lambda, Pxhbarxhinv)
                do j = 1, nvar+1
                    ! First term: tau * ws[H[j]] (element-wise, not matrix multiplication!)
                    obj_row(j) = tau * ws(H_subsample(j))
                    ! Second term: sum_i lambda[i] * Pxhbarxhinv[i, j]
                    do i = 1, n_hbar
                        obj_row(j) = obj_row(j) + lambda(i) * Pxhbarxhinv(i, j)
                    end do
                end do

                ! Store objective row at ms+1 (last row of tableau)
                do j = 1, nvar+1
                    gammaxs(ms + 1, j) = obj_row(j)
                end do


            end if

            ! Set bs for objective row at ms+1 (last row of tableau)
            bs(ms + 1) = 0.0d0

            ! Run PPRO-specific simplex on subsample
            ! NOTE: Pass ms (subsample size WITH aggregates), not n_subsample

            ! FIX: Copy to correctly-sized array to avoid column-major dimension mismatch
            do i = 1, ms+1
                do j = 1, nvar+1
                    gammaxs_simplex(i, j) = gammaxs(i, j)
                end do
                bs_simplex(i) = bs(i)
            end do

            ! Pass m+1-sized work arrays because the simplex stores the objective row
            ! at row mv+1 even when mv == m.
            remaining = maxit - iter_total
            if (remaining <= 0) then
                call set_failure(2)
                return
            end if
            call run_simplex_ppro(gammaxs_simplex, bs_simplex, IBs, freevarrows, r1s, r2s, rr, ws, &
                                  m + 1, ms, nvar, tau, tol, remaining, bland, iter_attempt, no_pivot_flag, &
                                  simplex_converged)
            iter_total = iter_total + iter_attempt
            if (no_pivot_flag .or. (.not. simplex_converged)) then
                if (force_full_sample) then
                    call set_failure(2)
                    return
                end if
                empty_pivot_count = empty_pivot_count + 1
                mmm = mmm * 2.0d0
                not_new_sl_sh = .true.
                if (empty_pivot_count >= max_empty_pivot_retries) then
                    force_full_sample = .true.
                end if
                cycle
            end if


            ! Extract solution from subsample
            ! NOTE: Use ms for extraction too
            ! Use bs_simplex which has the correct values after simplex
            call extract_solution(gammaxs_simplex, bs_simplex, IBs, ms, nvar, estimate, &
                                  u_subsample, v_subsample, r1s)


            ll_candidate = estimate(1) + estimate(2) * z_sorted(rd)
            d_ll_candidate = estimate(2)

            ! Compute raw full residuals for candidate verification.
            do i = 1, m
                r_raw(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
            end do

            call map_candidate_H(h_map_ok, h_failure_code)
            if (.not. h_map_ok) then
                call set_failure(h_failure_code)
                return
            end if

            if (always_same_h_refit) then
                call try_same_h_refit(refit_recertified)
                if (.not. refit_recertified) then
                    do i = 1, m
                        r_raw(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
                    end do
                end if
            end if


            ! BUG FIX: Do NOT sort H_mat! R code does not sort it.
            ! The order returned by the simplex (r1 indices) must be preserved
            ! for warm-start to work correctly.
            ! if (H_mat_sorted(rd, 1) > H_mat_sorted(rd, 2)) then
            !     j = H_mat_sorted(rd, 1)
            !     H_mat_sorted(rd, 1) = H_mat_sorted(rd, 2)
            !     H_mat_sorted(rd, 2) = j
            ! end if


            ! Check bad signs and certify using raw residuals.
            if (count(sl) > 0 .or. count(sh) > 0) then
                n_bad_signs = 0
                do i = 1, m
                    if ((sh(i) .and. r_raw(i) <= 0.0d0) .or. (sl(i) .and. r_raw(i) >= 0.0d0)) then
                        n_bad_signs = n_bad_signs + 1
                    end if
                end do

                if (n_bad_signs > 0) then
                    ! R logic: if bad.signs > 0.1 * ms, double mmm; otherwise remove bad obs
                    if (dble(n_bad_signs) > 0.1d0 * dble(ms)) then
                        ! Too many bad signs: double M and retry
                        mmm = mmm * 2.0d0
                        not_new_sl_sh = .true.
                        ! Continue while loop - will rebuild with larger M
                    else
                        ! Few bad signs: remove them from sl/sh and retry
                        do i = 1, m
                            if (sh(i) .and. r_raw(i) <= 0.0d0) sh(i) = .false.
                            if (sl(i) .and. r_raw(i) >= 0.0d0) sl(i) = .false.
                        end do
                        not_new_sl_sh = .false.
                        ! Continue while loop - will rebuild with adjusted sl/sh
                    end if
                else
                    accept_subsample = certify_llqr_candidate(H_candidate, r_raw, A, m, nvar, res_tol)
                    if ((.not. accept_subsample) .and. (.not. always_same_h_refit)) then
                        call try_same_h_refit(refit_recertified)
                        if (refit_recertified .and. n_bad_signs == 0) accept_subsample = .true.
                    end if
                    if (n_bad_signs > 0) then
                        if (dble(n_bad_signs) > 0.1d0 * dble(ms)) then
                            mmm = mmm * 2.0d0
                            not_new_sl_sh = .true.
                        else
                            do i = 1, m
                                if (sh(i) .and. r_raw(i) <= 0.0d0) sh(i) = .false.
                                if (sl(i) .and. r_raw(i) >= 0.0d0) sl(i) = .false.
                            end do
                            not_new_sl_sh = .false.
                        end if
                    else if (accept_subsample) then
                        ll_est_sorted(rd) = ll_candidate
                        d_ll_est_sorted(rd) = d_ll_candidate
                        H_mat_sorted(rd, :) = H_candidate
                        r = r_raw
                        do i = 1, nvar+1
                            r(H_candidate(i)) = 0.0d0
                        end do
                        residual_prev = r
                        call store_current_cache()
                        not_optimal = .false.
                    else
                        call handle_certification_reject(cert_reject_return)
                        if (cert_reject_return) return
                    end if
                end if
            else
                accept_subsample = certify_llqr_candidate(H_candidate, r_raw, A, m, nvar, res_tol)
                if ((.not. accept_subsample) .and. (.not. always_same_h_refit)) then
                    call try_same_h_refit(refit_recertified)
                    if (refit_recertified .and. n_bad_signs == 0) accept_subsample = .true.
                end if
                if (accept_subsample) then
                    ll_est_sorted(rd) = ll_candidate
                    d_ll_est_sorted(rd) = d_ll_candidate
                    H_mat_sorted(rd, :) = H_candidate
                    r = r_raw
                    do i = 1, nvar+1
                        r(H_candidate(i)) = 0.0d0
                    end do
                    residual_prev = r
                    call store_current_cache()
                    not_optimal = .false.
                else
                    call handle_certification_reject(cert_reject_return)
                    if (cert_reject_return) return
                end if
            end if

        else if (rd > 2) then
            ! ========================================================
            ! Step 3 (rd>2): Incremental warm start (R lines 664-732)
            ! ========================================================
            ! Reuse xhinv from previous round (R line 666)
            ! R: xhinv <- gammaxs.temp[1:(nvar + 1),]
            do i = 1, nvar+1
                do j = 1, nvar+1
                    xhinv(i, j) = xhinv_stored(i, j)
                end do
            end do

            ! Map H from previous round to subsample coordinates (same as rd==2)
            do i = 1, nvar+1
                H_subsample(i) = 0
                do j = 1, n_subsample
                    if (idx_not_jl_or_jh(j) == H_prev(i)) then
                        H_subsample(i) = j
                        exit
                    end if
                end do
                if (H_subsample(i) == 0) then
                    use_independent_init = .true.
                    independent_trigger = 2
                    not_new_sl_sh = .true.
                    cycle attempt_loop
                end if
            end do


            ! Identify positive and negative residuals in subsample (same as rd==2)
            n_idpos = 0
            n_idneg = 0
            do i = 1, n_subsample
                ! Check if this index is in H
                is_in_H = .false.
                do j = 1, nvar+1
                    if (H_subsample(j) == i) then
                        is_in_H = .true.
                        exit
                    end if
                end do

                if (.not. is_in_H) then
                    r_idx = idx_not_jl_or_jh(i)
                    if (residual_prev(r_idx) > 0.0d0) then
                        n_idpos = n_idpos + 1
                        idpos(n_idpos) = i
                    else if (residual_prev(r_idx) < 0.0d0) then
                        n_idneg = n_idneg + 1
                        idneg(n_idneg) = i
                    end if
                end if
            end do

            ! Build Hbar (same as rd==2)
            n_hbar = n_idpos + n_idneg
            do i = 1, n_idpos
                Hbar(i) = idpos(i)
            end do
            do i = 1, n_idneg
                Hbar(n_idpos + i) = idneg(i)
            end do
            if (has_sl_agg .and. has_sh_agg) then
                Hbar(n_hbar + 1) = m + 1
                Hbar(n_hbar + 2) = m + 2
                n_hbar = n_hbar + 2
            else if (has_sl_agg) then
                Hbar(n_hbar + 1) = m + 1
                n_hbar = n_hbar + 1
            else if (has_sh_agg) then
                Hbar(n_hbar + 1) = m + 2
                n_hbar = n_hbar + 1
            end if

            ! Reuse stored basis block exactly (top nvar+1 rows)
            ! R: gammaxs.temp[1:(nvar+1),] <- gammaxs.temp[1:(nvar+1),]
            !    bs.temp[1:(nvar+1)] <- bs.temp[1:(nvar+1)]
            ! NOTE: For rd>2, these are the previous LP's basis rows - DO NOT recompute!
            do i = 1, nvar+1
                do j = 1, nvar+1
                    gammaxs_temp(i, j) = xhinv_stored(i, j)
                end do
                bs_temp(i) = bs_stored(i)
            end do


            ! ========================================================
            ! Incremental update of positive Hbar rows (R lines 668-682)
            ! ========================================================
            ! Get current positive residual indices in original data
            ! R: idx_Hbar_pos2 <- idx_not_jl_or_jh[idpos]
            do i = 1, n_idpos
                idx_Hbar_pos2(i) = idx_not_jl_or_jh(idpos(i))
            end do

            ! Check which ones match previous round
            n_matched_pos = 0
            n_unmatched_pos = 0
            do i = 1, n_idpos
                ! R: matched_rows_pos <- idx_Hbar_pos2 %in% idx_Hbar_pos
                matched_rows_pos(i) = .false.
                do j = 1, n_pos_prev
                    if (idx_Hbar_pos2(i) == idx_Hbar_pos(j)) then
                        matched_rows_pos(i) = .true.
                        row_map_pos(i) = j  ! Store which row in gammaxs_pos to reuse
                        n_matched_pos = n_matched_pos + 1
                        exit
                    end if
                end do
                if (.not. matched_rows_pos(i)) then
                    n_unmatched_pos = n_unmatched_pos + 1
                    unmatched_idx_pos(n_unmatched_pos) = i
                end if
            end do


            ! Copy matched rows from previous gammaxs_pos
            ! R lines 673-675: indices <- (nvar+2):(nvar+1+length(u.in.IBs))
            !                  gammaxs.temp[indices[rows.pos2], ] <- gammaxs.pos[rows.pos, ]
            do i = 1, n_idpos
                curr_idx = nvar + 1 + i  ! Position in gammaxs_temp
                if (matched_rows_pos(i)) then
                    ! Reuse stored row
                    do j = 1, nvar+1
                        gammaxs_temp(curr_idx, j) = gammaxs_pos(row_map_pos(i), j)
                    end do
                    bs_temp(curr_idx) = bs_pos(row_map_pos(i))
                else
                    ! Compute new row: -A[idx,:] %*% xhinv (R lines 679-681)
                    ! R: Pxhbarxhinv.pos <- A[rows.pos2.nm,] %*% xhinv
                    ! R: gammaxs.temp[indices[!matched_rows_pos],] <- - Pxhbarxhinv.pos
                    r_idx = idx_Hbar_pos2(i)
                    do j = 1, nvar+1
                        gammaxs_temp(curr_idx, j) = 0.0d0
                        do k = 1, nvar+1
                            gammaxs_temp(curr_idx, j) = gammaxs_temp(curr_idx, j) - A(r_idx, k) * xhinv(k, j)
                        end do
                    end do
                    ! R: bs.temp[indices[!matched_rows_pos]] <- - Pxhbarxhinv.pos %*% y[idx_not_jl_or_jh[H]] + y[rows.pos2.nm]
                    ! Since gammaxs_temp = -Pxhbarxhinv.pos, we want gammaxs_temp %*% y[H] + y[rows]
                    bs_temp(curr_idx) = 0.0d0
                    do j = 1, nvar+1
                        bs_temp(curr_idx) = bs_temp(curr_idx) + gammaxs_temp(curr_idx, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                    end do
                    bs_temp(curr_idx) = bs_temp(curr_idx) + y(r_idx)

                end if
            end do

            ! ========================================================
            ! Incremental update of negative Hbar rows (R lines 684-698)
            ! ========================================================
            do i = 1, n_idneg
                idx_Hbar_neg2(i) = idx_not_jl_or_jh(idneg(i))
            end do

            n_matched_neg = 0
            n_unmatched_neg = 0
            do i = 1, n_idneg
                matched_rows_neg(i) = .false.
                do j = 1, n_neg_prev
                    if (idx_Hbar_neg2(i) == idx_Hbar_neg(j)) then
                        matched_rows_neg(i) = .true.
                        row_map_neg(i) = j
                        n_matched_neg = n_matched_neg + 1
                        exit
                    end if
                end do
                if (.not. matched_rows_neg(i)) then
                    n_unmatched_neg = n_unmatched_neg + 1
                    unmatched_idx_neg(n_unmatched_neg) = i
                end if
            end do

            ! Copy matched rows or compute new ones
            do i = 1, n_idneg
                curr_idx = nvar + 1 + n_idpos + i  ! Position in gammaxs_temp
                if (matched_rows_neg(i)) then
                    do j = 1, nvar+1
                        gammaxs_temp(curr_idx, j) = gammaxs_neg(row_map_neg(i), j)
                    end do
                    bs_temp(curr_idx) = bs_neg(row_map_neg(i))
                else
                    ! R: Pxhbarxhinv.neg <- - A[rows.neg2.nm,] %*% xhinv (line 695)
                    r_idx = idx_Hbar_neg2(i)
                    do j = 1, nvar+1
                        gammaxs_temp(curr_idx, j) = 0.0d0
                        do k = 1, nvar+1
                            gammaxs_temp(curr_idx, j) = gammaxs_temp(curr_idx, j) + A(r_idx, k) * xhinv(k, j)
                        end do
                    end do
                    ! R line 697: bs.temp[...] <- - Pxhbarxhinv.neg %*% y[...] - y[rows.neg2.nm]
                    ! Pxhbarxhinv.neg = -A*xhinv, gammaxs_temp = A*xhinv
                    ! R wants: -Pxhbarxhinv.neg*y[H] - y[rows] = A*xhinv*y[H] - y[rows] = gammaxs_temp*y[H] - y[rows]
                    bs_temp(curr_idx) = 0.0d0
                    do j = 1, nvar+1
                        bs_temp(curr_idx) = bs_temp(curr_idx) + gammaxs_temp(curr_idx, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                    end do
                    bs_temp(curr_idx) = bs_temp(curr_idx) - y(r_idx)

                end if
            end do

            ! ========================================================
            ! Update aggregate rows (R lines 700-724)
            ! ========================================================
            if (has_sl_agg .and. has_sh_agg) then
                ! R lines 701-704: Update both sl and sh aggregates
                ! gammaxs.temp[m+1,] <- gammaxs.temp[m+1,] %*% xhinv
                do j = 1, nvar+1
                    temp_vec(j) = 0.0d0
                    do k = 1, nvar+1
                        temp_vec(j) = temp_vec(j) + gammaxs_temp(m+1, k) * xhinv(k, j)
                    end do
                end do
                do j = 1, nvar+1
                    gammaxs_temp(m+1, j) = temp_vec(j)
                end do
                ! bs.temp[m+1] <- gammaxs.temp[m+1,] %*% y[idx_not_jl_or_jh[H]] - bs.temp[m+1]
                min_k = 0.0d0
                do j = 1, nvar+1
                    min_k = min_k + gammaxs_temp(m+1, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                end do
                bs_temp(m+1) = min_k - bs_temp(m+1)

                ! gammaxs.temp[m+2,] <- - gammaxs.temp[m+2,] %*% xhinv
                do j = 1, nvar+1
                    temp_vec(j) = 0.0d0
                    do k = 1, nvar+1
                        temp_vec(j) = temp_vec(j) - gammaxs_temp(m+2, k) * xhinv(k, j)
                    end do
                end do
                do j = 1, nvar+1
                    gammaxs_temp(m+2, j) = temp_vec(j)
                end do
                ! bs.temp[m+2] <- gammaxs.temp[m+2,] %*% y[...] + bs.temp[m+2]
                min_k = 0.0d0
                do j = 1, nvar+1
                    min_k = min_k + gammaxs_temp(m+2, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                end do
                bs_temp(m+2) = min_k + bs_temp(m+2)

            else if (has_sl_agg) then
                ! Only sl aggregate (R lines 709-710)
                do j = 1, nvar+1
                    temp_vec(j) = 0.0d0
                    do k = 1, nvar+1
                        temp_vec(j) = temp_vec(j) + gammaxs_temp(m+1, k) * xhinv(k, j)
                    end do
                end do
                do j = 1, nvar+1
                    gammaxs_temp(m+1, j) = temp_vec(j)
                end do
                min_k = 0.0d0
                do j = 1, nvar+1
                    min_k = min_k + gammaxs_temp(m+1, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                end do
                bs_temp(m+1) = min_k - bs_temp(m+1)

            else if (has_sh_agg) then
                ! Only sh aggregate (R lines 715-716)
                do j = 1, nvar+1
                    temp_vec(j) = 0.0d0
                    do k = 1, nvar+1
                        temp_vec(j) = temp_vec(j) - gammaxs_temp(m+2, k) * xhinv(k, j)
                    end do
                end do
                do j = 1, nvar+1
                    gammaxs_temp(m+2, j) = temp_vec(j)
                end do
                min_k = 0.0d0
                do j = 1, nvar+1
                    min_k = min_k + gammaxs_temp(m+2, j) * y(idx_not_jl_or_jh(H_subsample(j)))
                end do
                bs_temp(m+2) = min_k + bs_temp(m+2)

            end if

            ! ========================================================
            ! Build final gammaxs and bs from gammaxs_temp (R lines 705-724)
            ! ========================================================
            ! For rd>2, gammaxs_temp contains:
            !   Rows 1:(nvar+1): xhinv
            !   Rows (nvar+2):(nvar+1+n_idpos): positive Hbar rows
            !   Rows (nvar+2+n_idpos):(nvar+1+n_idpos+n_idneg): negative Hbar rows
            !   Rows m+1, m+2: aggregates (if they exist)
            ! R: ms <- nvar + 1 + length(u.in.IBs) + length(v.in.IBs) + (1 if sl) + (1 if sh)
            ! R: gammaxs <- gammaxs.temp[c(1:(ms-2),m+1,m+2),]
            n_xhinv_hbar = nvar + 1 + n_idpos + n_idneg  ! This is ms-2 when both aggregates exist


            if (has_sl_agg .and. has_sh_agg) then
                ! Copy xhinv and Hbar rows 1:(ms-2), then aggregates m+1, m+2
                do i = 1, n_xhinv_hbar
                    do j = 1, nvar+1
                        gammaxs(i, j) = gammaxs_temp(i, j)
                    end do
                    bs(i) = bs_temp(i)
                end do
                do j = 1, nvar+1
                    gammaxs(n_xhinv_hbar + 1, j) = gammaxs_temp(m+1, j)
                    gammaxs(n_xhinv_hbar + 2, j) = gammaxs_temp(m+2, j)
                end do
                bs(n_xhinv_hbar + 1) = bs_temp(m+1)
                bs(n_xhinv_hbar + 2) = bs_temp(m+2)
                ! Build lambda (R line 707)
                do i = 1, n_idpos
                    lambda(i) = tau * ws(idpos(i))
                end do
                do i = 1, n_idneg
                    lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                end do
                lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
                lambda(n_idpos + n_idneg + 2) = tau
            else if (has_sl_agg) then
                ! Copy xhinv and Hbar rows 1:(ms-1), then aggregate m+1
                do i = 1, n_xhinv_hbar
                    do j = 1, nvar+1
                        gammaxs(i, j) = gammaxs_temp(i, j)
                    end do
                    bs(i) = bs_temp(i)
                end do
                do j = 1, nvar+1
                    gammaxs(n_xhinv_hbar + 1, j) = gammaxs_temp(m+1, j)
                end do
                bs(n_xhinv_hbar + 1) = bs_temp(m+1)
                do i = 1, n_idpos
                    lambda(i) = tau * ws(idpos(i))
                end do
                do i = 1, n_idneg
                    lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                end do
                lambda(n_idpos + n_idneg + 1) = 1.0d0 - tau
            else if (has_sh_agg) then
                ! Copy xhinv and Hbar rows 1:(ms-1), then aggregate m+2
                do i = 1, n_xhinv_hbar
                    do j = 1, nvar+1
                        gammaxs(i, j) = gammaxs_temp(i, j)
                    end do
                    bs(i) = bs_temp(i)
                end do
                do j = 1, nvar+1
                    gammaxs(n_xhinv_hbar + 1, j) = gammaxs_temp(m+2, j)
                end do
                bs(n_xhinv_hbar + 1) = bs_temp(m+2)
                do i = 1, n_idpos
                    lambda(i) = tau * ws(idpos(i))
                end do
                do i = 1, n_idneg
                    lambda(n_idpos + i) = (1.0d0 - tau) * ws(idneg(i))
                end do
                lambda(n_idpos + n_idneg + 1) = tau
            else
                ! No aggregates - copy all xhinv and Hbar rows (1:ms)
                do i = 1, n_xhinv_hbar
                    do j = 1, nvar+1
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

            ! Compute Pxhbarxhinv and objective row (R lines 726-731)
            ! R: Pxhbarxhinv <- - gammaxs[(nvar+2):ms,]
            ! This extracts the Hbar rows (which are already stored in gammaxs after xhinv rows)
            ! R: gammaxs <- rbind(gammaxs, tau * ws[H] + t(lambda) %*% Pxhbarxhinv)
            ! Compute current ms (number of rows in gammaxs after assembly)
            if (has_sl_agg .and. has_sh_agg) then
                ms_current = n_xhinv_hbar + 2
                n_lambda = n_idpos + n_idneg + 2
            else if (has_sl_agg .or. has_sh_agg) then
                ms_current = n_xhinv_hbar + 1
                n_lambda = n_idpos + n_idneg + 1
            else
                ms_current = n_xhinv_hbar
                n_lambda = n_idpos + n_idneg
            end if

            ! Compute objective row
            do j = 1, nvar+1
                obj_row(j) = tau * ws(H_subsample(j))
                ! Pxhbarxhinv are rows (nvar+2):ms_current
                do i = 1, n_lambda
                    ! Pxhbarxhinv is -gammaxs[nvar+1+i,]
                    obj_row(j) = obj_row(j) - lambda(i) * gammaxs(nvar + 1 + i, j)
                end do
            end do
            do j = 1, nvar+1
                gammaxs(ms_current + 1, j) = obj_row(j)
            end do
            bs(ms_current + 1) = 0.0d0


            ! Initialize IBs, freevarrows, r1s, r2s (same as rd==2)
            do i = 1, nvar+1
                IBs(i) = i
            end do
            do i = 1, n_idpos
                u_in_IBs(i) = idpos(i) + nvar + 1
                IBs(nvar + 1 + i) = u_in_IBs(i)
            end do
            do i = 1, n_idneg
                v_in_IBs(i) = idneg(i) + nvar + 1 + ms
                IBs(nvar + 1 + n_idpos + i) = v_in_IBs(i)
            end do
            if (has_sl_agg .and. has_sh_agg) then
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + 2*ms - 1
                IBs(nvar + 1 + n_idpos + n_idneg + 2) = nvar + 1 + ms
            else if (has_sl_agg) then
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + 2*ms
            else if (has_sh_agg) then
                IBs(nvar + 1 + n_idpos + n_idneg + 1) = nvar + 1 + ms
            end if

            n_hbar_core = n_idpos + n_idneg
            do i = 1, ms + 1
                freevarrows(i) = .false.
            end do
            do i = 1, nvar+1
                freevarrows(i) = .true.
            end do
            do i = nvar + 2, nvar + 1 + n_hbar_core
                freevarrows(i) = .false.
            end do
            if (has_sl_agg .and. has_sh_agg) then
                freevarrows(nvar + 1 + n_hbar_core + 1) = .true.  ! sl aggregate row
                freevarrows(nvar + 1 + n_hbar_core + 2) = .true.  ! sh aggregate row
            else if (has_sl_agg .or. has_sh_agg) then
                freevarrows(nvar + 1 + n_hbar_core + 1) = .true.  ! single aggregate row
            end if
            freevarrows(ms + 1) = .true.

            do i = 1, nvar+1
                r1s(i) = H_subsample(i) + nvar + 1
                r2s(i) = r1s(i) + ms
            end do

            ! Copy to simplex arrays (same as rd==2)
            do i = 1, ms+1
                do j = 1, nvar+1
                    gammaxs_simplex(i, j) = gammaxs(i, j)
                end do
                bs_simplex(i) = bs(i)
            end do

            ! Run simplex (same as rd==2)
            remaining = maxit - iter_total
            if (remaining <= 0) then
                call set_failure(2)
                return
            end if
            call run_simplex_ppro(gammaxs_simplex, bs_simplex, IBs, freevarrows, r1s, r2s, rr, ws, &
                                  m + 1, ms, nvar, tau, tol, remaining, bland, iter_attempt, no_pivot_flag, &
                                  simplex_converged)
            iter_total = iter_total + iter_attempt
            if (no_pivot_flag .or. (.not. simplex_converged)) then
                if (force_full_sample) then
                    call set_failure(2)
                    return
                end if
                empty_pivot_count = empty_pivot_count + 1
                mmm = mmm * 2.0d0
                not_new_sl_sh = .true.
                if (empty_pivot_count >= max_empty_pivot_retries) then
                    force_full_sample = .true.
                end if
                cycle
            end if

            ! Extract solution (same as rd==2)
            call extract_solution(gammaxs_simplex, bs_simplex, IBs, ms, nvar, estimate, &
                                  u_subsample, v_subsample, r1s)

            ll_candidate = estimate(1) + estimate(2) * z_sorted(rd)
            d_ll_candidate = estimate(2)

            ! Compute raw full residuals for candidate verification.
            do i = 1, m
                r_raw(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
            end do

            call map_candidate_H(h_map_ok, h_failure_code)
            if (.not. h_map_ok) then
                call set_failure(h_failure_code)
                return
            end if

            if (always_same_h_refit) then
                call try_same_h_refit(refit_recertified)
                if (.not. refit_recertified) then
                    do i = 1, m
                        r_raw(i) = y(i) - (A(i, 1) * estimate(1) + A(i, 2) * estimate(2))
                    end do
                end if
            end if

            ! Check bad signs and certify using raw residuals.
            if (count(sl) > 0 .or. count(sh) > 0) then
                n_bad_signs = 0
                do i = 1, m
                    if ((sh(i) .and. r_raw(i) <= 0.0d0) .or. (sl(i) .and. r_raw(i) >= 0.0d0)) then
                        n_bad_signs = n_bad_signs + 1
                    end if
                end do

                if (n_bad_signs > 0) then
                    if (dble(n_bad_signs) > 0.1d0 * dble(ms)) then
                        mmm = mmm * 2.0d0
                        not_new_sl_sh = .true.
                    else
                        do i = 1, m
                            if (sh(i) .and. r_raw(i) <= 0.0d0) sh(i) = .false.
                            if (sl(i) .and. r_raw(i) >= 0.0d0) sl(i) = .false.
                        end do
                        not_new_sl_sh = .false.
                    end if
                else
                    accept_subsample = certify_llqr_candidate(H_candidate, r_raw, A, m, nvar, res_tol)
                    if ((.not. accept_subsample) .and. (.not. always_same_h_refit)) then
                        call try_same_h_refit(refit_recertified)
                        if (refit_recertified .and. n_bad_signs == 0) accept_subsample = .true.
                    end if
                    if (n_bad_signs > 0) then
                        if (dble(n_bad_signs) > 0.1d0 * dble(ms)) then
                            mmm = mmm * 2.0d0
                            not_new_sl_sh = .true.
                        else
                            do i = 1, m
                                if (sh(i) .and. r_raw(i) <= 0.0d0) sh(i) = .false.
                                if (sl(i) .and. r_raw(i) >= 0.0d0) sl(i) = .false.
                            end do
                            not_new_sl_sh = .false.
                        end if
                    else if (accept_subsample) then
                        ll_est_sorted(rd) = ll_candidate
                        d_ll_est_sorted(rd) = d_ll_candidate
                        H_mat_sorted(rd, :) = H_candidate
                        r = r_raw
                        do i = 1, nvar+1
                            r(H_candidate(i)) = 0.0d0
                        end do
                        residual_prev = r
                        call store_current_cache()
                        not_optimal = .false.
                    else
                        call handle_certification_reject(cert_reject_return)
                        if (cert_reject_return) return
                    end if
                end if
            else
                accept_subsample = certify_llqr_candidate(H_candidate, r_raw, A, m, nvar, res_tol)
                if ((.not. accept_subsample) .and. (.not. always_same_h_refit)) then
                    call try_same_h_refit(refit_recertified)
                    if (refit_recertified .and. n_bad_signs == 0) accept_subsample = .true.
                end if
                if (accept_subsample) then
                    ll_est_sorted(rd) = ll_candidate
                    d_ll_est_sorted(rd) = d_ll_candidate
                    H_mat_sorted(rd, :) = H_candidate
                    r = r_raw
                    do i = 1, nvar+1
                        r(H_candidate(i)) = 0.0d0
                    end do
                    residual_prev = r
                    call store_current_cache()
                    not_optimal = .false.
                else
                    call handle_certification_reject(cert_reject_return)
                    if (cert_reject_return) return
                end if
            end if

        end if  ! rd == 2 or rd > 2

        end do attempt_loop  ! while (not_optimal)

    end do  ! rd = 2, rounds

    ! ============================================================
    ! Unsort results back to original order (matching R track_order=TRUE)
    ! R code: if (track_order) ll_est <- ll_est[order(original_order)]
    ! ============================================================
    do rd = 1, rounds
        k = z_order(rd)  ! Original position of this sorted element
        ll_est(k) = ll_est_sorted(rd)
        d_ll_est(k) = d_ll_est_sorted(rd)
        do i = 1, nvar+1
            H_mat(k, i) = H_mat_sorted(rd, i)
        end do
    end do

contains

    double precision function llqr_kernel_weight(u_val, case_val, pi_val)
        implicit none
        double precision, intent(in) :: u_val, pi_val
        integer, intent(in) :: case_val

        if (case_val == 2) then
            if (abs(u_val) <= 1.0d0) then
                llqr_kernel_weight = 0.75d0 * (1.0d0 - u_val * u_val)
            else
                llqr_kernel_weight = 0.0d0
            end if
        else
            llqr_kernel_weight = exp(-0.5d0 * u_val * u_val) / sqrt(2.0d0 * pi_val)
        end if
    end function llqr_kernel_weight

    subroutine validate_previous_H(H_idx, valid)
        implicit none
        integer, intent(in) :: H_idx(nvar+1)
        logical, intent(out) :: valid
        integer :: vi, vj

        valid = .true.
        do vi = 1, nvar + 1
            if (H_idx(vi) < 1 .or. H_idx(vi) > m) then
                valid = .false.
                return
            end if
            do vj = vi + 1, nvar + 1
                if (H_idx(vi) == H_idx(vj)) then
                    valid = .false.
                    return
                end if
            end do
        end do
    end subroutine validate_previous_H

    subroutine reorder_cold_tableau(theta_offset, n_individual, n_total, success)
        implicit none
        double precision, intent(in) :: theta_offset(nvar+1)
        integer, intent(in) :: n_individual, n_total
        logical, intent(out) :: success
        double precision :: gx_tmp(m+1, nvar+1), bv_tmp(m+1)
        integer :: IB_tmp(m+1)
        logical :: fvr_tmp(m+1), used_row(m)
        integer :: ci, cj, src, dest

        success = .false.
        used_row = .false.

        ! Put coefficient-basic rows first and convert their shifted values
        ! back to absolute coefficients. Aggregate rows are protected, so a
        ! coefficient-basic row must be an individual row.
        do ci = 1, nvar + 1
            src = 0
            do cj = 1, n_total
                if (IBs(cj) == ci) then
                    src = cj
                    exit
                end if
            end do
            if (src < 1 .or. src > n_individual) return
            used_row(src) = .true.
            do cj = 1, nvar + 1
                gx_tmp(ci, cj) = gammaxs_simplex(src, cj)
            end do
            bv_tmp(ci) = bs_simplex(src) + theta_offset(ci)
            IB_tmp(ci) = IBs(src)
            fvr_tmp(ci) = freevarrows(src)
        end do

        dest = nvar + 1
        do src = 1, n_individual
            if (.not. used_row(src)) then
                dest = dest + 1
                do cj = 1, nvar + 1
                    gx_tmp(dest, cj) = gammaxs_simplex(src, cj)
                end do
                bv_tmp(dest) = bs_simplex(src)
                IB_tmp(dest) = IBs(src)
                fvr_tmp(dest) = freevarrows(src)
            end if
        end do
        do src = n_individual + 1, n_total
            dest = dest + 1
            do cj = 1, nvar + 1
                gx_tmp(dest, cj) = gammaxs_simplex(src, cj)
            end do
            bv_tmp(dest) = bs_simplex(src)
            IB_tmp(dest) = IBs(src)
            fvr_tmp(dest) = freevarrows(src)
        end do
        if (dest /= n_total) return

        do cj = 1, nvar + 1
            gx_tmp(n_total + 1, cj) = gammaxs_simplex(n_total + 1, cj)
        end do
        bv_tmp(n_total + 1) = bs_simplex(n_total + 1)
        IB_tmp(n_total + 1) = IBs(n_total + 1)
        fvr_tmp(n_total + 1) = freevarrows(n_total + 1)

        do ci = 1, n_total + 1
            do cj = 1, nvar + 1
                gammaxs_simplex(ci, cj) = gx_tmp(ci, cj)
            end do
            bs_simplex(ci) = bv_tmp(ci)
            IBs(ci) = IB_tmp(ci)
            freevarrows(ci) = fvr_tmp(ci)
        end do
        success = .true.
    end subroutine reorder_cold_tableau

    subroutine run_shifted_reduced_initialization(success)
        implicit none
        logical, intent(out) :: success
        double precision :: shifted_rhs, delta_est(nvar+1)
        integer :: si, sj, src, total_rows
        logical :: reorder_success

        success = .false.
        total_rows = ms

        do si = 1, total_rows
            if (si <= n_subsample) then
                src = si
            else if (has_sl_agg .and. si == n_subsample + 1) then
                src = m + 1
            else
                src = m + 2
            end if

            shifted_rhs = bs_temp(src)
            do sj = 1, nvar + 1
                shifted_rhs = shifted_rhs - gammaxs_temp(src, sj) * theta_prev(sj)
            end do

            if (shifted_rhs < 0.0d0) then
                do sj = 1, nvar + 1
                    gammaxs_simplex(si, sj) = -gammaxs_temp(src, sj)
                end do
                bs_simplex(si) = -shifted_rhs
                IBs(si) = nvar + 1 + total_rows + si
            else
                do sj = 1, nvar + 1
                    gammaxs_simplex(si, sj) = gammaxs_temp(src, sj)
                end do
                bs_simplex(si) = shifted_rhs
                IBs(si) = nvar + 1 + si
            end if
            freevarrows(si) = (si > n_subsample)
        end do

        IBs(total_rows + 1) = 0
        freevarrows(total_rows + 1) = .true.
        bs_simplex(total_rows + 1) = 0.0d0
        do sj = 1, nvar + 1
            gammaxs_simplex(total_rows + 1, sj) = 0.0d0
            do si = 1, total_rows
                if (IBs(si) > nvar + 1 .and. IBs(si) <= nvar + 1 + total_rows) then
                    gammaxs_simplex(total_rows + 1, sj) = gammaxs_simplex(total_rows + 1, sj) - &
                        tau * ws(si) * gammaxs_simplex(si, sj)
                else if (IBs(si) > nvar + 1 + total_rows) then
                    gammaxs_simplex(total_rows + 1, sj) = gammaxs_simplex(total_rows + 1, sj) - &
                        (1.0d0 - tau) * ws(si) * gammaxs_simplex(si, sj)
                end if
            end do
        end do

        do si = 1, nvar + 1
            r1s(si) = si
            r2s(si) = 0
        end do
        rr = 0.0d0

        remaining = maxit - iter_total
        if (remaining <= 0) return
        call run_simplex_full_llqr(gammaxs_simplex, bs_simplex, IBs, freevarrows, r1s, r2s, rr, ws, &
                                   m + 1, total_rows, nvar, tau, tol, remaining, bland, iter_attempt, &
                                   no_pivot_flag, simplex_converged)
        iter_total = iter_total + iter_attempt
        if (no_pivot_flag .or. (.not. simplex_converged)) return

        call extract_solution(gammaxs_simplex, bs_simplex, IBs, total_rows, nvar, delta_est, &
                              u_subsample, v_subsample, r1s)
        do si = 1, nvar + 1
            estimate(si) = theta_prev(si) + delta_est(si)
        end do

        call map_candidate_H(h_map_ok, h_failure_code)
        if (.not. h_map_ok) return

        call reorder_cold_tableau(theta_prev, n_subsample, total_rows, reorder_success)
        if (.not. reorder_success) return
        success = .true.
    end subroutine run_shifted_reduced_initialization

    subroutine run_full_active_recovery(success)
        implicit none
        logical, intent(out) :: success
        double precision :: zero_offset(nvar+1)
        integer :: fi, fj
        logical :: reorder_success

        success = .false.
        n_subsample = m
        ms = m
        has_sl_agg = .false.
        has_sh_agg = .false.
        sl = .false.
        sh = .false.
        do fi = 1, m
            idx_not_jl_or_jh(fi) = fi
            ws(fi) = w(fi)
            if (y(fi) < 0.0d0) then
                do fj = 1, nvar + 1
                    gammaxs_simplex(fi, fj) = -A(fi, fj)
                end do
                bs_simplex(fi) = -y(fi)
                IBs(fi) = nvar + 1 + m + fi
            else
                do fj = 1, nvar + 1
                    gammaxs_simplex(fi, fj) = A(fi, fj)
                end do
                bs_simplex(fi) = y(fi)
                IBs(fi) = nvar + 1 + fi
            end if
            freevarrows(fi) = .false.
        end do

        IBs(m + 1) = 0
        freevarrows(m + 1) = .true.
        bs_simplex(m + 1) = 0.0d0
        do fj = 1, nvar + 1
            gammaxs_simplex(m + 1, fj) = 0.0d0
            do fi = 1, m
                if (IBs(fi) > nvar + 1 .and. IBs(fi) <= nvar + 1 + m) then
                    gammaxs_simplex(m + 1, fj) = gammaxs_simplex(m + 1, fj) - &
                        tau * w(fi) * gammaxs_simplex(fi, fj)
                else
                    gammaxs_simplex(m + 1, fj) = gammaxs_simplex(m + 1, fj) - &
                        (1.0d0 - tau) * w(fi) * gammaxs_simplex(fi, fj)
                end if
            end do
        end do

        do fi = 1, nvar + 1
            r1s(fi) = fi
            r2s(fi) = 0
        end do
        rr = 0.0d0

        remaining = maxit - iter_total
        if (remaining <= 0) return
        call run_simplex_full_llqr(gammaxs_simplex, bs_simplex, IBs, freevarrows, r1s, r2s, rr, w, &
                                   m + 1, m, nvar, tau, tol, remaining, bland, iter_attempt, &
                                   no_pivot_flag, simplex_converged)
        iter_total = iter_total + iter_attempt
        if (no_pivot_flag .or. (.not. simplex_converged)) return

        call extract_solution(gammaxs_simplex, bs_simplex, IBs, m, nvar, estimate, &
                              u_subsample, v_subsample, r1s)
        call map_candidate_H(h_map_ok, h_failure_code)
        if (.not. h_map_ok) return

        zero_offset = 0.0d0
        call reorder_cold_tableau(zero_offset, m, m, reorder_success)
        if (.not. reorder_success) return
        success = .true.
    end subroutine run_full_active_recovery

    logical function terminal_cert_failure()
        implicit none

        terminal_cert_failure = (ms >= m) .or. ((.not. has_sl_agg) .and. (.not. has_sh_agg)) .or. force_full_sample
    end function terminal_cert_failure

    subroutine handle_certification_reject(should_return)
        implicit none
        logical, intent(out) :: should_return

        if (terminal_cert_failure()) then
            call set_failure(1)
            should_return = .true.
        else
            mmm = mmm * 2.0d0
            not_new_sl_sh = .true.
            should_return = .false.
        end if
    end subroutine handle_certification_reject

    subroutine store_current_cache()
        implicit none
        integer :: si, sj, sms_org, scurr_idx

        do si = 1, nvar+1
            do sj = 1, nvar+1
                xhinv_stored(si, sj) = gammaxs_simplex(si, sj)
            end do
            bs_stored(si) = bs_simplex(si)
        end do

        sms_org = ms
        if (has_sl_agg) sms_org = sms_org - 1
        if (has_sh_agg) sms_org = sms_org - 1

        n_pos_prev = 0
        n_neg_prev = 0
        do si = nvar+2, sms_org
            if (IBs(si) > nvar + 1 + ms) then
                n_neg_prev = n_neg_prev + 1
                scurr_idx = IBs(si) - nvar - 1 - ms
                idx_Hbar_neg(n_neg_prev) = idx_not_jl_or_jh(scurr_idx)
                do sj = 1, nvar+1
                    gammaxs_neg(n_neg_prev, sj) = gammaxs_simplex(si, sj)
                end do
                bs_neg(n_neg_prev) = bs_simplex(si)
            else
                n_pos_prev = n_pos_prev + 1
                scurr_idx = IBs(si) - nvar - 1
                idx_Hbar_pos(n_pos_prev) = idx_not_jl_or_jh(scurr_idx)
                do sj = 1, nvar+1
                    gammaxs_pos(n_pos_prev, sj) = gammaxs_simplex(si, sj)
                end do
                bs_pos(n_pos_prev) = bs_simplex(si)
            end if
        end do
    end subroutine store_current_cache

    subroutine set_failure(code)
        implicit none
        integer, intent(in) :: code

        ierr = code
        if (rd >= 1 .and. rd <= rounds) then
            failed_eval = rd
            H_mat_sorted(rd, :) = 0
        else
            failed_eval = 1
        end if
    end subroutine set_failure

    subroutine map_candidate_H(success, failure_code)
        implicit none
        logical, intent(out) :: success
        integer, intent(out) :: failure_code
        integer :: mi, mj, reduced_idx

        success = .false.
        failure_code = 0
        H_candidate = 0

        do mi = 1, nvar+1
            reduced_idx = r1s(mi) - nvar - 1
            if (reduced_idx < 1 .or. reduced_idx > n_subsample) then
                failure_code = 5
                return
            end if
            H_candidate(mi) = idx_not_jl_or_jh(reduced_idx)
            if (H_candidate(mi) < 1 .or. H_candidate(mi) > m) then
                failure_code = 3
                return
            end if
        end do

        do mi = 1, nvar+1
            do mj = mi + 1, nvar+1
                if (H_candidate(mi) == H_candidate(mj)) then
                    failure_code = 3
                    return
                end if
            end do
        end do

        success = .true.
    end subroutine map_candidate_H

    logical function certify_llqr_candidate(H_idx, r_vec, A_mat, m_loc, nvar_loc, res_tol_loc)
        implicit none
        integer, intent(in) :: m_loc, nvar_loc
        integer, intent(in) :: H_idx(nvar_loc+1)
        double precision, intent(in) :: r_vec(m_loc), A_mat(m_loc, nvar_loc+1), res_tol_loc
        integer :: i, j
        double precision :: det2

        certify_llqr_candidate = .true.
        do i = 1, nvar_loc + 1
            if (H_idx(i) < 1 .or. H_idx(i) > m_loc) then
                certify_llqr_candidate = .false.
                return
            end if
            do j = i + 1, nvar_loc + 1
                if (H_idx(i) == H_idx(j)) then
                    certify_llqr_candidate = .false.
                    return
                end if
            end do
        end do

        if (nvar_loc == 1) then
            det2 = A_mat(H_idx(1), 1) * A_mat(H_idx(2), 2) - A_mat(H_idx(1), 2) * A_mat(H_idx(2), 1)
            if (abs(det2) <= 1.0d-10) then
                certify_llqr_candidate = .false.
                return
            end if
        end if

        do i = 1, nvar_loc + 1
            if (abs(r_vec(H_idx(i))) > res_tol_loc) then
                certify_llqr_candidate = .false.
                return
            end if
        end do
    end function certify_llqr_candidate

    subroutine try_same_h_refit(recertified)
        implicit none
        logical, intent(out) :: recertified
        double precision :: xh_refit(nvar+1, nvar+1)
        double precision :: xhinv_refit(nvar+1, nvar+1)
        double precision :: estimate_refit(nvar+1)
        double precision :: det2_refit
        logical :: inv_success_refit
        integer :: ri, rj

        recertified = .false.

        do ri = 1, nvar + 1
            if (H_candidate(ri) < 1 .or. H_candidate(ri) > m) return
            do rj = ri + 1, nvar + 1
                if (H_candidate(ri) == H_candidate(rj)) return
            end do
        end do

        if (nvar /= 1) return

        do ri = 1, nvar + 1
            do rj = 1, nvar + 1
                xh_refit(ri, rj) = A(H_candidate(ri), rj)
            end do
        end do

        det2_refit = xh_refit(1, 1) * xh_refit(2, 2) - xh_refit(1, 2) * xh_refit(2, 1)
        if (abs(det2_refit) <= 1.0d-10) return

        call inv22(xh_refit, xhinv_refit, inv_success_refit)
        if (.not. inv_success_refit) return

        do ri = 1, nvar + 1
            estimate_refit(ri) = 0.0d0
            do rj = 1, nvar + 1
                estimate_refit(ri) = estimate_refit(ri) + xhinv_refit(ri, rj) * y(H_candidate(rj))
            end do
        end do

        do ri = 1, m
            r_raw(ri) = y(ri) - (A(ri, 1) * estimate_refit(1) + A(ri, 2) * estimate_refit(2))
        end do

        n_bad_signs = 0
        do ri = 1, m
            if ((sh(ri) .and. r_raw(ri) <= 0.0d0) .or. (sl(ri) .and. r_raw(ri) >= 0.0d0)) then
                n_bad_signs = n_bad_signs + 1
            end if
        end do

        if (.not. certify_llqr_candidate(H_candidate, r_raw, A, m, nvar, res_tol)) return

        estimate = estimate_refit
        do ri = 1, nvar + 1
            bs_simplex(ri) = estimate_refit(ri)
        end do
        ll_candidate = estimate(1) + estimate(2) * z_sorted(rd)
        d_ll_candidate = estimate(2)
        recertified = .true.
    end subroutine try_same_h_refit

    subroutine run_simplex_full_llqr(gx, bv, IBv, fvr, r1v, r2v, rrv, wv, ldgx, mv, nvr, &
                                     tv, tl, mxit, bld, iters, no_pivot, converged)
        implicit none
        integer, intent(in) :: ldgx, mv, nvr, mxit
        double precision, intent(inout) :: gx(ldgx, nvr+1), bv(mv+1)
        integer, intent(inout) :: IBv(mv+1), r1v(nvr+1), r2v(nvr+1)
        logical, intent(inout) :: fvr(mv+1)
        double precision, intent(inout) :: rrv(2, nvr+1)
        double precision, intent(in) :: wv(mv), tv, tl
        logical, intent(in) :: bld
        integer, intent(out) :: iters
        logical, intent(out) :: no_pivot, converged

        double precision :: yyv(mv+1), eev(mv+1), k_valsv(mv+1)
        integer :: ii, jj, kk, t_rrv, tsepv, tv_val
        integer :: idx_offset
        double precision :: rrlv, min_kv, pivot_val

        iters = 0
        no_pivot = .false.
        converged = .false.
        t_rrv = 1
        tsepv = 1
        tv_val = 0

        do while (iters < mxit)
            do ii = 1, nvr+1
                rrv(1, ii) = gx(mv+1, ii)
                if (r2v(ii) /= 0) then
                    idx_offset = r1v(ii) - 1 - nvr
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

            if (bld) then
                if (any(rrv(1,:) < -tl)) then
                    tv_val = huge(1)
                    do ii = 1, nvr+1
                        if (rrv(1,ii) < -tl .and. r1v(ii) < tv_val) then
                            tv_val = r1v(ii)
                            t_rrv = ii
                            tsepv = 1
                        end if
                    end do
                else
                    tv_val = huge(1)
                    do ii = 1, nvr+1
                        if (rrv(2,ii) < -tl .and. r2v(ii) < tv_val) then
                            tv_val = r2v(ii)
                            t_rrv = ii
                            tsepv = 2
                        end if
                    end do
                end if
            else
                do jj = 1, nvr+1
                    do ii = 1, 2
                        if (abs(rrv(ii,jj) - rrlv) < tl) then
                            t_rrv = jj
                            tsepv = ii
                            if (tsepv == 1) then
                                tv_val = r1v(t_rrv)
                            else
                                tv_val = r2v(t_rrv)
                            end if
                            goto 200
                        end if
                    end do
                end do
200             continue
            end if

            if (r2v(t_rrv) /= 0) then
                if (tsepv == 1) then
                    do ii = 1, mv+1
                        yyv(ii) = gx(ii, t_rrv)
                    end do
                else
                    do ii = 1, mv+1
                        yyv(ii) = -gx(ii, t_rrv)
                    end do
                end if

                min_kv = huge(1.0d0)
                kk = 0
                do ii = 1, mv+1
                    if (yyv(ii) > tl .and. .not. fvr(ii)) then
                        k_valsv(ii) = bv(ii) / yyv(ii)
                        if (k_valsv(ii) < min_kv - tl) then
                            min_kv = k_valsv(ii)
                            kk = ii
                        else if (abs(k_valsv(ii) - min_kv) < tl .and. bld) then
                            if (kk == 0 .or. IBv(ii) < IBv(kk)) kk = ii
                        end if
                    end if
                end do

                if (kk == 0) then
                    no_pivot = .true.
                    exit
                end if

                if (tsepv /= 1) then
                    idx_offset = r1v(t_rrv) - 1 - nvr
                    if (idx_offset >= 1 .and. idx_offset <= mv) then
                        yyv(mv+1) = yyv(mv+1) + wv(idx_offset)
                    end if
                end if
            else
                do ii = 1, mv+1
                    yyv(ii) = gx(ii, t_rrv)
                end do

                min_kv = huge(1.0d0)
                kk = 0
                if (yyv(mv+1) < 0.0d0) then
                    do ii = 1, mv+1
                        if (yyv(ii) > tl .and. .not. fvr(ii)) then
                            k_valsv(ii) = bv(ii) / yyv(ii)
                            if (k_valsv(ii) < min_kv - tl) then
                                min_kv = k_valsv(ii)
                                kk = ii
                            else if (abs(k_valsv(ii) - min_kv) < tl .and. bld) then
                                if (kk == 0 .or. IBv(ii) < IBv(kk)) kk = ii
                            end if
                        end if
                    end do
                else
                    do ii = 1, mv+1
                        if (yyv(ii) < -tl .and. .not. fvr(ii)) then
                            k_valsv(ii) = -bv(ii) / yyv(ii)
                            if (k_valsv(ii) < min_kv - tl) then
                                min_kv = k_valsv(ii)
                                kk = ii
                            else if (abs(k_valsv(ii) - min_kv) < tl .and. bld) then
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

            if (IBv(kk) <= (mv + nvr + 1)) then
                do ii = 1, mv+1
                    gx(ii, t_rrv) = 0.0d0
                end do
                gx(kk, t_rrv) = 1.0d0
                r1v(t_rrv) = IBv(kk)
                r2v(t_rrv) = IBv(kk) + mv
            else
                do ii = 1, mv+1
                    gx(ii, t_rrv) = 0.0d0
                end do
                gx(kk, t_rrv) = -1.0d0
                idx_offset = IBv(kk) - mv - nvr - 1
                if (idx_offset >= 1 .and. idx_offset <= mv) then
                    gx(mv+1, t_rrv) = wv(idx_offset)
                end if
                r1v(t_rrv) = IBv(kk) - mv
                r2v(t_rrv) = IBv(kk)
            end if

            do jj = 1, nvr+1
                pivot_val = gx(kk, jj)
                do ii = 1, mv+1
                    gx(ii, jj) = gx(ii, jj) - eev(ii) * pivot_val
                end do
            end do

            pivot_val = bv(kk)
            do ii = 1, mv+1
                bv(ii) = bv(ii) - eev(ii) * pivot_val
            end do

            IBv(kk) = tv_val
            iters = iters + 1
        end do
    end subroutine run_simplex_full_llqr

    ! PPRO-specific simplex algorithm (matches R's llqr_tau_seq_ppro lines 731-831)
    subroutine run_simplex_ppro(gx, bv, IBv, fvr, r1v, r2v, rrv, wv, ldgx, mv, nvr, &
                                tv, tl, mxit, bld, iters, no_pivot, converged)
        implicit none
        integer, intent(in) :: ldgx, mv, nvr, mxit
        double precision, intent(inout) :: gx(ldgx, nvr+1), bv(mv+1)
        integer, intent(inout) :: IBv(mv+1), r1v(nvr+1), r2v(nvr+1)
        logical, intent(inout) :: fvr(mv+1)
        double precision, intent(inout) :: rrv(2, nvr+1)
        double precision, intent(in) :: wv(mv), tv, tl
        logical, intent(in) :: bld
        integer, intent(out) :: iters
        logical, intent(out) :: no_pivot, converged

        ! Local variables
        double precision :: yyv(mv+1), eev(mv+1), k_valsv(mv+1)
        integer :: ii, jj, kk, t_rrv, tsepv, tv_val
        integer :: best_var_id, best_col
        integer :: idx_offset
        double precision :: rrlv, min_kv, pivot_val

        iters = 0
        no_pivot = .false.
        converged = .false.
        t_rrv = 1
        tsepv = 1
        tv_val = 0


        do while (iters < mxit)
            ! Step 2: Compute reduced costs (PPRO version)
            ! R code line 735-736:
            ! rr[1, ] <- gammaxs[ms+1, ]
            ! rr[2, ] <- (ws[r1 - 1 - nvar] - rr[1, ])


            do ii = 1, nvr+1
                rrv(1, ii) = gx(mv+1, ii)

                ! PPRO formula: ws[r1 - 1 - nvar] - rr[1, ]
                idx_offset = r1v(ii) - 1 - nvr
                if (idx_offset >= 1 .and. idx_offset <= mv) then
                    rrv(2, ii) = wv(idx_offset) - rrv(1, ii)
                else
                    rrv(2, ii) = -rrv(1, ii)  ! If out of range, just negate
                end if

                ! R code: rr[1, r2==0] <- -abs(rr[1, r2==0])
                ! When r2[i] == 0, negate rr[1,i] to make it negative
                if (r2v(ii) == 0) then
                    rrv(1, ii) = -abs(rrv(1, ii))
                end if
            end do


            ! Check optimality
            rrlv = minval(rrv)
            if (rrlv >= -tl) then
                converged = .true.
                exit
            end if

            ! Step 3: Choose entering variable
            if (bld) then
                ! Bland's rule: match R by choosing the smallest eligible variable id.
                best_var_id = huge(1)
                best_col = 0
                if (any(rrv(1,:) < -tl)) then
                    do jj = 1, nvr+1
                        if (rrv(1,jj) < -tl .and. r1v(jj) < best_var_id) then
                            best_var_id = r1v(jj)
                            best_col = jj
                        end if
                    end do
                    tsepv = 1
                else
                    do jj = 1, nvr+1
                        if (rrv(2,jj) < -tl .and. r2v(jj) < best_var_id) then
                            best_var_id = r2v(jj)
                            best_col = jj
                        end if
                    end do
                    tsepv = 2
                end if
                t_rrv = best_col
                tv_val = best_var_id
            else
                ! Standard rule: match R's which(rr == min(rr), arr.ind=TRUE)[1,]
                ! scan order for a 2 x p matrix: column outer, row inner.
                do jj = 1, nvr+1
                    do ii = 1, 2
                        if (rrv(ii,jj) == rrlv) then
                            t_rrv = jj
                            tsepv = ii
                            if (tsepv == 1) then
                                tv_val = r1v(t_rrv)
                            else
                                tv_val = r2v(t_rrv)
                            end if
                            goto 100
                        end if
                    end do
                end do
100             continue
            end if


            ! Step 4 & 5: Choose leaving variable (PPRO version)
            ! R code line 776-780
            ! IMPORTANT: Only copy rows 1:mv+1, not all ldgx rows!
            if (tsepv == 1) then
                do ii = 1, mv+1
                    yyv(ii) = gx(ii, t_rrv)
                end do
            else
                do ii = 1, mv+1
                    yyv(ii) = -gx(ii, t_rrv)
                end do
            end if


            ! Ratio test (R code line 783-800)
            ! k <- bs / yy
            ! k <- which((min(k[yy > 0 & !freevarrow]) == k) & (yy > 0))
            min_kv = huge(1.0d0)
            kk = 0

            do ii = 1, mv+1
                if (yyv(ii) > tl .and. .not. fvr(ii)) then
                    k_valsv(ii) = bv(ii) / yyv(ii)
                    if (k_valsv(ii) < min_kv - tl) then
                        min_kv = k_valsv(ii)
                        kk = ii
                    else if (abs(k_valsv(ii) - min_kv) < tl .and. bld) then
                        ! Bland's rule tie-breaking: choose smallest index
                        if (IBv(ii) < IBv(kk)) then
                            kk = ii
                        end if
                    end if
                end if
            end do

            if (kk == 0) then
                no_pivot = .true.
                exit
            end if


            ! Step 6': Adjust yy if entering variable is v_i (R code line 805-807)
            ! if (tsep != 1){ yy[ms + 1] <- yy[ms + 1] + ws[r1[t_rr] - 1 - nvar] }
            if (tsepv /= 1) then
                idx_offset = r1v(t_rrv) - 1 - nvr
                if (idx_offset >= 1 .and. idx_offset <= mv) then
                    yyv(mv+1) = yyv(mv+1) + wv(idx_offset)
                end if
            end if

            ! Step 6: Pivoting (R code line 809-828)
            ! ee <- yy / yy[k]; ee[k] <- 1 - 1 / yy[k]
            do ii = 1, mv+1
                if (ii == kk) then
                    eev(ii) = 1.0d0 - 1.0d0 / yyv(kk)
                else
                    eev(ii) = yyv(ii) / yyv(kk)
                end if
            end do

            ! Update pivot column and r1, r2 (R code line 812-823)
            if (IBv(kk) <= (mv + nvr + 1)) then
                ! R: if (IBs[k] <= (ms + nvar + 1))
                ! IMPORTANT: Only zero out rows 1:mv+1, not all ldgx rows!
                do ii = 1, mv+1
                    gx(ii, t_rrv) = 0.0d0
                end do
                gx(kk, t_rrv) = 1.0d0
                r1v(t_rrv) = IBv(kk)
                r2v(t_rrv) = IBv(kk) + mv
            else
                ! R: else { gammaxs[k, t_rr] <- -1; gammaxs[(ms + 1), t_rr] <- ws[IBs[k] - ms - nvar - 1] }
                ! IMPORTANT: Only zero out rows 1:mv+1, not all ldgx rows!
                do ii = 1, mv+1
                    gx(ii, t_rrv) = 0.0d0
                end do
                gx(kk, t_rrv) = -1.0d0
                idx_offset = IBv(kk) - mv - nvr - 1
                if (idx_offset >= 1 .and. idx_offset <= mv) then
                    gx(mv+1, t_rrv) = wv(idx_offset)
                end if
                r1v(t_rrv) = IBv(kk) - mv
                r2v(t_rrv) = IBv(kk)
            end if

            ! Update all columns (tcrossprod)
            do jj = 1, nvr+1
                pivot_val = gx(kk, jj)
                do ii = 1, mv+1
                    gx(ii, jj) = gx(ii, jj) - eev(ii) * pivot_val
                end do
            end do

            ! Update b
            pivot_val = bv(kk)
            do ii = 1, mv+1
                bv(ii) = bv(ii) - eev(ii) * pivot_val
            end do

            IBv(kk) = tv_val

            ! R code: if (t <= nvar+1) { freevarrow[k] <- TRUE }
            ! Mark row as free variable if a coefficient entered
            if (tv_val <= nvr + 1) then
                fvr(kk) = .true.
            end if

            iters = iters + 1

        end do

    end subroutine run_simplex_ppro

    ! Extract solution from basis
    subroutine extract_solution(gx, bv, IBv, mv, nvr, est, uv, vv, r1v)
        implicit none
        integer, intent(in) :: mv, nvr
        double precision, intent(in) :: gx(mv+1, nvr+1), bv(mv+1)
        integer, intent(in) :: IBv(mv+1), r1v(nvr+1)
        double precision, intent(out) :: est(nvr+1), uv(mv), vv(mv)
        integer :: ii, jj, u_idx, v_idx
        logical :: u_tmp(mv), v_tmp(mv)

        ! Initialize
        est = 0.0d0
        uv = 0.0d0
        vv = 0.0d0

        ! Extract from basis (only loop through m, not m+1)
        do ii = 1, mv
            if (IBv(ii) > nvr + 1 .and. IBv(ii) <= nvr + 1 + mv) then
                ! u variable is basic
                u_idx = IBv(ii) - nvr - 1
                if (u_idx >= 1 .and. u_idx <= mv) then
                    uv(u_idx) = bv(ii)
                end if
            else if (IBv(ii) > nvr + 1 + mv) then
                ! v variable is basic
                v_idx = IBv(ii) - nvr - 1 - mv
                if (v_idx >= 1 .and. v_idx <= mv) then
                    vv(v_idx) = bv(ii)
                end if
            else if (IBv(ii) >= 1 .and. IBv(ii) <= nvr + 1) then
                ! Coefficient variable is basic
                est(IBv(ii)) = bv(ii)
            end if
        end do
    end subroutine extract_solution

end subroutine llqr_ppro_fortran
