module mod_time_integ
    implicit none
    
contains
    subroutine rk_explicit(order,dt,u,N)
        implicit none
        integer(4), intent(in)  :: order
        real(8), intent(in)     :: dt
        real(8)                 :: c_i(4),a_ij(4),b_i(4)
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
        ! starting the loop over RK
        do i = 1,order
            ! compute arguments
            beta = t_n + c_i(i) * dt
            alpha(:) = u(:) + dt * a_ij(i) * g_i(:)
            ! Computing g_i
            g_i(:) = R(alpha,beta)
            ! Accummulate F = sum(b_i*g_i)
            F(:) = F(:) + b_i(i) * g_i(:)
        end do
        ! Update the solution to the next time-step
        u(:) = u(:) + dt*F(:)
    end subroutine rk_explicit
end module mod_time_integ