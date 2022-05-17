module mod_time_integ
    implicit none
    
contains
    subroutine rk_explicit(order,c1,Re,u,N)
        implicit none
        ! integer(4), parameter :: N = 20, tmax = 2950 ! N=20 case
        ! ! integer(4), parameter :: N = 100, tmax = 3200 ! N=100 case
        ! integer(4) :: Re, k, p, q, t
        ! real(8) :: c1, dt, total_time, real_part, imaginary_part, E(N), start, finish
        ! complex(8) :: onei, up_iq_uq, uk(N), ut(N), err(tmax), sum_error(tmax), Ri_0(N), ui_1(N), Ri_1(N)
        integer(4), intent(in)    :: order, N, Re ! N is a parameter?
        real(8), intent(in)       :: c1
        integer(4)                :: i,k
        real(8)                   :: dt, c_i(4),a_ij(4),b_i(4), real_part, imaginary_part
        complex(8)                :: beta(N), alpha(N), g_i(N), F(N)
        complex(8), intent(inout) :: u(N)
        dt = c1 * Re / N ** 2.0
        ! create RK table based on order
        if (order == 1) then
            c_i = 0.0
            a_ij = 0.0
            b_i = 0.0
            b_i(1) = 1.0
        else if (order == 2) then
            c_i = 0.0
            a_ij = 0.0
            b_i = 0.0
            c_i(2) = 0.0
            a_ij(2) = 0.5
            b_i(1) = 0.5
            b_i(2) = 0.5
        else if (order == 3) then
            print *, "Invalid input, no output available"
            error stop
        else if (order == 4) then
            c_i = 0.0
            a_ij = 0.0
            b_i = 0.0
            c_i(2) = 0.5
            c_i(3) = 0.5
            c_i(4) = 1.0
            a_ij(2) = 0.5
            a_ij(3) = 0.5
            a_ij(4) = 1.0
            b_i(1) = 1.0/6.0
            b_i(2) = 1.0/3.0
            b_i(3) = 1.0/3.0
            b_i(4) = 1.0/6.0
        else
            print *, "Invalid input, no output available"
            error stop
        end if
        ! Initialization
        F(:) = 0.0
        g_i(:) = 0.0
        do k = 1,N
            real_part = 1.0/k ! k ** (-1.0)
            imaginary_part = 0.0
            u(k) = cmplx(real_part, imaginary_part, 8)
        end do
        ! starting the loop over RK
        do i = 1,order
            ! compute arguments
            beta(:) = t_n + c_i(i) * dt
            alpha(:) = u(:) + dt * a_ij(i) * g_i(:)
            ! Computing g_i
            g_i(:) = R(beta, alpha)
            ! Accummulate F = sum(b_i*g_i)
            F(:) = F(:) + b_i(i) * g_i(:)
        end do
        ! Update the solution to the next time-step
        u(:) = u(:) + dt*F(:)
    end subroutine rk_explicit
end module mod_time_integ