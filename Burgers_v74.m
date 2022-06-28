% This program solves the Burgers equation using Fourier analysis
% and explicit Runge-Kutta 4th order method
%            -------
% written by Tariq Ridwan: https://tariqridwan.github.io/
% Barcelona Supercomputing Center // Universitat Politècnica de Catalunya
%% Burgers version 7-4: RK4 steady-state calculation
close all
clear
%% Physical description
Re = 40; % Reynolds number
N = 20; % total number of fourier modes
c1 = 0.473;  % courant coeff
tmax = 620; % total time
% N = 100; % total number of fourier modes
% c1 = 1.288;  % courant coeff
% tmax = 1250; % total time
dt = c1 * Re / N^2;
total_time = dt*tmax;
%% Initial condition
uk = zeros(1,N);
for k = 1:N   % uk(1) to uk(20)
    uk(k) = 1/k + 0i;
end
%% solving uk for different modes
Ri_0 = zeros(1,N); % for RK 1st stage
ui_1 = zeros(1,N);
Ri_1 = zeros(1,N); % for RK 2nd stage
ui_2 = zeros(1,N);
Ri_2 = zeros(1,N); % for RK 3rd stage
ui_3 = zeros(1,N);
Ri_3 = zeros(1,N); % for RK 4th stage
ut = zeros(1,N);
err = zeros(1,tmax);
ut(1) = uk(1);      % ut(1) = uk(1) as u(1) is constant
ui_1(1) = uk(1);      % ui_1(1) = uk(1) as u(1) is constant
ui_2(1) = uk(1);      % ui_1(1) = uk(1) as u(1) is constant
ui_3(1) = uk(1);      % ui_1(1) = uk(1) as u(1) is constant
tic
for t = 1:tmax
    for k = 2:N   % ut(2) to ut(20)
        up_iq_uq = 0i; % RK 1st stage starts
        for p = -N:N
            q = k - p;
            if q >= -N && q <= N
                if p >= 1 && q >= 1
                    up_iq_uq = up_iq_uq + uk(p)*(0+1i)*q*uk(q);
                elseif q < 0
                    up_iq_uq = up_iq_uq + uk(p)*(0+1i)*q*conj(uk(-q));
                elseif p < 0
                    up_iq_uq = up_iq_uq + conj(uk(-p))*(0+1i)*q*uk(q);
                end
            end
        end
%         ut(k) = uk(k) + dt*( -up_iq_uq - k^2/Re*uk(k) ); % Euler normal
        Ri_0(k) = -up_iq_uq - k^2/Re*uk(k);
        ui_1(k) = uk(k) + dt/2*Ri_0(k); % RK 1st stage ends
    end
    
    for k = 2:N   % ut(2) to ut(20)        
        up_iq_uq = 0i; % RK 2nd stage starts
        for p = -N:N
            q = k - p;
            if q >= -N && q <= N
                if p >= 1 && q >= 1
                    up_iq_uq = up_iq_uq + ui_1(p)*(0+1i)*q*ui_1(q);
                elseif q < 0
                    up_iq_uq = up_iq_uq + ui_1(p)*(0+1i)*q*conj(ui_1(-q));
                elseif p < 0
                    up_iq_uq = up_iq_uq + conj(ui_1(-p))*(0+1i)*q*ui_1(q);
                end
            end
        end
        Ri_1(k) = -up_iq_uq - k^2/Re*ui_1(k);
        ui_2(k) = uk(k) + dt/2*Ri_1(k); % RK 2nd stage ends
    end
    
    for k = 2:N   % ut(2) to ut(20)
        up_iq_uq = 0i; % RK 3rd stage starts
        for p = -N:N
            q = k - p;
            if q >= -N && q <= N
                if p >= 1 && q >= 1
                    up_iq_uq = up_iq_uq + ui_2(p)*(0+1i)*q*ui_2(q);
                elseif q < 0
                    up_iq_uq = up_iq_uq + ui_2(p)*(0+1i)*q*conj(ui_2(-q));
                elseif p < 0
                    up_iq_uq = up_iq_uq + conj(ui_2(-p))*(0+1i)*q*ui_2(q);
                end
            end
        end
        Ri_2(k) = -up_iq_uq - k^2/Re*ui_2(k);
        ui_3(k) = uk(k) + dt*Ri_2(k); % RK 3rd stage ends
    end
    
    for k = 2:N   % ut(2) to ut(20)
        up_iq_uq = 0i; % RK 4th stage starts
        for p = -N:N
            q = k - p;
            if q >= -N && q <= N
                if p >= 1 && q >= 1
                    up_iq_uq = up_iq_uq + ui_3(p)*(0+1i)*q*ui_3(q);
                elseif q < 0
                    up_iq_uq = up_iq_uq + ui_3(p)*(0+1i)*q*conj(ui_3(-q));
                elseif p < 0
                    up_iq_uq = up_iq_uq + conj(ui_3(-p))*(0+1i)*q*ui_3(q);
                end
            end
        end
        Ri_3(k) = -up_iq_uq - k^2/Re*ui_3(k); % RK 4th stage ends
    end

    for k = 2:N
%         ut(k) = uk(k) + dt*( -up_iq_uq - k^2/Re*uk(k) ); % Euler
%         ut(k) = uk(k) + dt*( Ri_0(k) + 0*Ri_1(k) + 0*Ri_2(k) + 0*Ri_3(k) ); % Euler normal
%         ut(k) = uk(k) + dt/2*( Ri_0(k) + 1*Ri_1(k) + 0*Ri_2(k) + 0*Ri_3(k) ); % RK2
        ut(k) = uk(k) + dt/6*( Ri_0(k) + 2*Ri_1(k) + 2*Ri_2(k) + Ri_3(k) ); % RK4
        err(t) = err(t) + ( ut(k) - uk(k) )^2;
        uk(k) = ut(k);
    end
end
toc
CPUtime = toc
%% Calculating total error
sum_error = sqrt(err); % to see Error VS time plot
%%  Plot Error VS time-step graph
figure(1)
% loglog(real(sum_error),'k.-')
semilogy(1:tmax,real(sum_error))
xlim([0 tmax])
ylim([1e-15 1e+0])
xlabel('time-step','interpreter','latex','fontsize',12)
ylabel('Error','interpreter','latex','fontsize',12)
grid on
% legend({'Error with increasing time-step'},'Interpreter','latex','fontsize',10,'Location','Northeast')
legend_plot = sprintf('RK4, $N=20$, $CFL$ = %.3f',c1);
% legend_plot = sprintf('RK4, $N=100$, $CFL$ = %.3f',c1);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Burgers equation: Error with increasing time-step','Interpreter','latex','fontsize',12)
% title_plot = sprintf('$N$ = %d',N);
% title(title_plot,'Interpreter','latex','fontsize',12)
print Burgers_ERROR-vs-Time-step.png -dpng
print Burgers_ERROR-vs-Time-step.eps -depsc
print Burgers_ERROR-vs-Time-step.pdf -dpdf
%%  Plot Error VS time graph
figure(2)
time_matrix = (dt:dt:total_time);
semilogy(time_matrix,real(sum_error))
xlim([0 time_matrix(tmax)])
ylim([1e-15 1e+0])
xlabel('time','interpreter','latex','fontsize',12)
ylabel('Error','interpreter','latex','fontsize',12)
grid on
% legend({'Error with increasing time'},'Interpreter','latex','fontsize',10,'Location','Northeast')
legend_plot = sprintf('RK4, $N=20$, $dt$ = %.4f',dt);
% legend_plot = sprintf('RK4, $N=100$, $dt$ = %.4f',dt);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Burgers equation: Error with increasing time','Interpreter','latex','fontsize',12)
% title_plot = sprintf('$N$ = %d',N);
% title(title_plot,'Interpreter','latex','fontsize',12)
print Burgers_ERROR-vs-Time.png -dpng
print Burgers_ERROR-vs-Time.eps -depsc
print Burgers_ERROR-vs-Time.pdf -dpdf
%% solving Ek for different modes
E = zeros(1,N);
for k = 1:N
    E(k) = ut(k)*conj(ut(k));
end
%% visualize Ek for N = 100
% figure(3)
% loglog(E,'ko-','LineWidth',1)
% xlim([1 1*N])
% ylim([1e-6 1e+0])
% xlabel('$k$','interpreter','latex','fontsize',12) 
% ylabel('$E_k$','interpreter','latex','fontsize',12)
% grid on
% % legend({'N=100'},'Interpreter','latex','fontsize',10,'Location','Northeast')
% legend_plot = sprintf('RK4, $N=100$, $CFL$ = %.3f',c1);
% legend(legend_plot,'Interpreter','latex','fontsize',11)
% title('Energy spectrum of Burgers equation''s steady-state solution','Interpreter','latex','fontsize',12)
% % title_plot = sprintf('$CFL$ = %.3f',c1);
% % title(title_plot,'Interpreter','latex','fontsize',12)
% print Burgers.png -dpng
% print Burgers.eps -depsc
% print Burgers.pdf -dpdf
% hold on
%% visualize Ek for N = 20
figure(3)
loglog(E,'+r-','LineWidth',1)
xlim([1 5*N])
ylim([1e-6 1e+0])
xlabel('$k$','interpreter','latex','fontsize',12) 
ylabel('$E_k$','interpreter','latex','fontsize',12)
grid on
% legend({'$N=20$'},'Interpreter','latex','fontsize',10,'Location','Northeast')
legend_plot = sprintf('RK4, $N=20$, $CFL$ = %.3f',c1);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Energy spectrum of Burgers equation''s steady-state solution','Interpreter','latex','fontsize',12)
% title_plot = sprintf('$CFL$ = %.3f',c1);
% title(title_plot,'Interpreter','latex','fontsize',12)
print Burgers.png -dpng
print Burgers.eps -depsc
print Burgers.pdf -dpdf
% hold on
%% Calculating & Plotting the real function 'u' over the real space
a_k = 2*real(ut);
b_k = 2*imag(ut);
u = zeros(1,10);
p = zeros(1,10);
% u_old = zeros(1,10);
for x = 1:151
u_old = 0;
    for k = 1:N
        p(x) = x-1;
        u(x) = u_old + a_k(k)*cos(k*(x-1)) + b_k(k)*sin(k*(x-1));
        u_old = u(x);
    end
end
figure(4)
plot(p,u,'.-')
% xlim([1 11])
ylim([-3 3])
% legend({'$u$'},'Interpreter','latex','fontsize',11,'Location','Best')
legend_plot = sprintf('$u$, when $N=%d$',N);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Steady-state solution of Burgers equation','Interpreter','latex','fontsize',12)
xlabel('$x$','interpreter','latex','fontsize',12) 
ylabel('$u$','interpreter','latex','fontsize',12)
grid on
print Burgers_u.png -dpng
print Burgers_u.eps -depsc
print Burgers_u.pdf -dpdf
savefig('Burgers_u.fig')