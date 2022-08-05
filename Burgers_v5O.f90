program Burgers_v5O
implicit none
! This program solves the Burgers equation using Fourier analysis
! and Euler explicit method
!            -------
! written by Tariq Ridwan: https://tariqridwan.github.io/
! Barcelona Supercomputing Center // Universitat Politècnica de Catalunya
! Burgers version 5.O: steady-state calculation_Optimized
!            -------
! Type declarations of different properties
    integer(4), parameter :: N = 20, tmax = 900 ! N=20 case
    ! integer(4), parameter :: N = 100, tmax = 1750 ! N=100 case
    integer(4) :: Re, k, p, q, t
    real(8) :: c1, dt, total_time, real_part, imaginary_part, E(N), start, finish
    complex(8) :: onei, up_iq_uq, uk(N), ut(N), err(tmax), sum_error(tmax)

! Give values to different parameters
    ! %% Physical description
    Re = 40 ! Reynolds number
    ! N = 20 ! total number of fourier modes
    ! tmax = 900 ! total time
    c1 = 0.3 ! courant coeff ! N=20 case
    ! c1 = 0.76 ! courant coeff ! N=100 case
    dt = c1 * Re / N ** 2.0
    total_time = dt*tmax
    k = N
    p = k
    q = k
    onei = cmplx(0.0, 1.0, 8)
    
! %% Initial condition
    ! uk = zeros(1,N); %% automatic in Fortran
    do k = 1,N
        real_part = 1.0/k ! k ** (-1.0)
        imaginary_part = 0.0
        uk(k) = cmplx(real_part, imaginary_part, 8) ! FLAGGGGGG :::: ????
        ! uk(k) = (real_part, imaginary_part)
        Print *, "INITIAL uk : ", k, uk(k)
    end do

! %% solving uk for different modes
    ut(1) = uk(1)  ! % ut(1) = uk(1) as u(1) is constant
    do k = 1,N
        print *, "INITIAL ut :  ", k, ut(k)
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
            ut(k) = uk(k) + dt*(-up_iq_uq-k**2.0/Re*uk(k))
            err(t) = err(t) + (ut(k)-uk(k))**2.0
            uk(k) = ut(k)
        end do
        ! print *, "t : ", t
    end do
    call cpu_time(finish) ! end calculating CPU-time 
    print '("Time = ",f6.4," seconds.")',finish-start ! print CPU-time
    ! print *, 'Time for second chunk of code was ', finish-start, 'seconds.'
    WRITE (*,*) 'Time of operation was ', finish-start, ' seconds'

    sum_error = sqrt(err)
    ! do t = 1,tmax
    !     Print *, "Error : ", sum_error(t) ! print all errors
    ! end do

    do k = 1,N
        Print *, "FINAL ut : ", k, ut(k) ! print all ut
    end do

! %% solving Ek for different modes
    do k = 1,N
        E(k) = ut(k)*conjg(ut(k))
        print *, "E : ", k, E(k) ! print all E
    end do

end program Burgers_v5O
