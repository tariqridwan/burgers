% This program solves the Burgers equation using Fourier analysis
% and explicit Runge-Kutta 2nd order method
%            -------
% written by Tariq Ridwan: tariq.ridwan@bsc.es
% Barcelona Supercomputing Center // Universitat Politècnica de Catalunya
%% Burgers version 7: RK2 steady-state calculation
close all
clear
%% Physical description
Re = 40; % Reynolds number
N = 20; % total number of fourier modes
c1 = 0.088;  % courant coeff
tmax = 2950; % total time
% N = 100; % total number of fourier modes
% c1 = 0.493;  % courant coeff
% tmax = 3200; % total time
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
ut = zeros(1,N);
err = zeros(1,tmax);
ut(1) = uk(1);      % ut(1) = uk(1) as u(1) is constant
ui_1(1) = uk(1);      % ui_1(1) = uk(1) as u(1) is constant
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
        up_iq_uq_1 = 0i; % RK 2nd stage starts
        for p = -N:N
            q = k - p;
            if q >= -N && q <= N
                if p >= 1 && q >= 1
                    up_iq_uq_1 = up_iq_uq_1 + ui_1(p)*(0+1i)*q*ui_1(q);
                elseif q < 0
                    up_iq_uq_1 = up_iq_uq_1 + ui_1(p)*(0+1i)*q*conj(ui_1(-q));
                elseif p < 0
                    up_iq_uq_1 = up_iq_uq_1 + conj(ui_1(-p))*(0+1i)*q*ui_1(q);
                end
            end
        end
        Ri_1(k) = -up_iq_uq_1 - k^2/Re*ui_1(k); % RK 2nd stage ends
%     end
% 
%     for k = 2:N
%         ut(k) = uk(k) + dt*( -up_iq_uq - k^2/Re*uk(k) ); % Euler normal
%         ut(k) = uk(k) + dt*( 1.0*Ri_0(k) + 0.0*Ri_1(k) ); % Euler
        ut(k) = uk(k) + dt*( 0.5*Ri_0(k) + 0.5*Ri_1(k) ); % RK2
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
% legend_plot = sprintf('RK2, $N=20$, $CFL$ = %.3f',c1);
legend_plot = sprintf('RK2, $N=100$, $CFL$ = %.3f',c1);
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
% legend_plot = sprintf('RK2, $N=20$, $dt$ = %.4f',dt);
legend_plot = sprintf('RK2, $N=100$, $dt$ = %.4f',dt);
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
% legend_plot = sprintf('RK2, $N=100$, $CFL$ = %.3f',c1);
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
legend_plot = sprintf('RK2, $N=20$, $CFL$ = %.3f',c1);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Energy spectrum of Burgers equation''s steady-state solution','Interpreter','latex','fontsize',12)
% title_plot = sprintf('$CFL$ = %.3f',c1);
% title(title_plot,'Interpreter','latex','fontsize',12)
print Burgers.png -dpng
print Burgers.eps -depsc
print Burgers.pdf -dpdf
hold on