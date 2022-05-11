program Burgers_v6O
implicit none
! Burgers version 5.O: steady-state calculation_Optimized

! Type declarations of different properties
    ! integer(4), parameter :: N = 20, tmax = 410 ! N=20 LES 4523
    integer(4), parameter :: N = 20, tmax = 1100 ! N=20 LES 05
    integer(4) :: Re, k, p, q, t, m
    real(8) :: c1, dt, total_time, real_part, imaginary_part, C_k, vt_star(N), E(N), nu_t(N), nu_eff(N), x, nut_inf, start, finish
    complex(8) :: onei, up_iq_uq, uk(N), ut(N), err(tmax), sum_error(tmax)
    ! nu_t(k) = nut_inf*sqrt( uk(N)*conjg(uk(N))/N )*vt_star(k)
    ! nu_eff(k) = 1.0/Re + nu_t(k) ! % for LES
    ! real(8) :: c1, dt, total_time, real_part, imaginary_part, C_k, x, start, finish
    ! complex(8) :: onei, up_iq_uq, uk(N), ut(N), err(tmax), sum_error(tmax), nu_t(N), nu_eff(N), nut_inf, vt_star(N), E(N)
    
! Give values to different parameters
    ! %% Physical description
    Re = 40 ! Reynolds number
    ! C_k = 0.4523
    ! c1 = 0.237 ! courant coeff ! N=20 LES 4523
    C_k = 0.05
    c1 = 0.056 ! courant coeff ! N=20 LES 05
    dt = c1 * Re / N ** 2.0
    total_time = dt*tmax
    k = N
    p = k
    q = k
    x = 1.0 ! to calculate e
    onei = cmplx(0.0, 1.0, 8)
    m = 2
    nut_inf = 0.31*(5-m)/(m+1)*sqrt(3.0-m)*C_k**(-3.0/2.0)
    print *, dt, total_time, nut_inf, C_k**(-3.0/2.0)
! %% Initial condition
    ! uk = zeros(1,N); %% automatic in Fortran
    do k = 1,N
        real_part = 1.0/k ! k ** (-1.0)
        imaginary_part = 0.0
        uk(k) = cmplx(real_part, imaginary_part, 8) ! FLAGGGGGG :::: ????
        ! uk(k) = (real_part, imaginary_part)
        Print *, "INITIAL uk : ", uk(k)
    end do

! %% solving uk for different modes
    ut(1) = uk(1)  ! % ut(1) = uk(1) as u(1) is constant
    do k = 1,N
        print *, "INITIAL ut :  ", ut(k)
    end do

    call cpu_time(start) ! start calculating CPU-time
    do t = 1,tmax
        do k = 2,N
            up_iq_uq = cmplx(0.0, 0.0, 8) ! (0.0, 0.0) %Initializing the term
            do p = -N,N ! % setting q values: -N < p < N
                q = k-p ! % Instead of q = -N:N
                if (q >= -N .and. q <= N) then
                    if (p >= 1 .and. q >= 1) then ! when p & q are positive
                        up_iq_uq = up_iq_uq + uk(p)*onei*q*uk(q)
                    else if (q < 0) then ! when only q is negative (p is positive)
                        up_iq_uq = up_iq_uq + uk(p)*onei*q*conjg(uk(-q))
                    else if (p < 0) then ! when only p is negative (q is positive)
                        up_iq_uq = up_iq_uq + conjg(uk(-p))*onei*q*uk(q)
                    else ! when p and or q = 0
                        up_iq_uq = up_iq_uq + 0.0
                    end if
                end if
            end do
            ! ut(k) = uk(k) + dt*(-up_iq_uq-k**2.0/Re*uk(k))
            vt_star(k) = 1.0 + 34.5*exp(x)**(-3.03*(N/k))
            nu_t(k) = nut_inf*sqrt( uk(N)*conjg(uk(N))/N )*vt_star(k)
            nu_eff(k) = 1.0/Re + nu_t(k) ! % for LES
            ! nu_eff(k) = 1/Re; ! % for normal & DNS
            ut(k) = uk(k) + dt*( -up_iq_uq - k**2.0 * nu_eff(k) * uk(k) )
            err(t) = err(t) + (ut(k)-uk(k))**2.0
            uk(k) = ut(k)
        end do
        ! print *, "t : ", t
    end do
    call cpu_time(finish) ! end calculating CPU-time
    ! print '("Time = ",f6.4," seconds.")',finish-start ! print CPU-time
    ! print '("Time = ",e5.6," seconds.")',finish-start ! print CPU-time
    print *, 'Time for second chunk of code was ', finish-start, 'seconds.'
    WRITE (*,*) 'Time of operation was ', finish-start, ' seconds'

    sum_error = sqrt(err)
    ! do t = 1,tmax
    !     Print *, "Error : ", sum_error(t) ! print all errors
    ! end do

    do k = 1,N
        Print *, "FINAL ut : ", ut(k) ! print all ut
    end do

! %% solving Ek for different modes
    do k = 1,N
        E(k) = ut(k)*conjg(ut(k))
        print *, "E : ", E(k) ! print all E
    end do

end program Burgers_v6O
