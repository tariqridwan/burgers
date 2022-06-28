% This program solves the Burgers equation using Fourier analysis
% and Euler explicit method
% with Large Eddy Simulation (LES) technique
%            -------
% written by Tariq Ridwan: https://tariqridwan.github.io/
% Barcelona Supercomputing Center // Universitat Politècnica de Catalunya
%% Burgers version 6.O: LES with steady-state calculation_Optimized
% format long
close all
clear
%% Physical description
Re = 40; % Reynolds number
N = 20; % total number of fourier modes
C_k = 0.4523;
c1 = 0.237;  % courant coeff
tmax = 410; % total time
% C_k = 0.305;
% c1 = 0.18;  % courant coeff
% tmax = 410; % total time
% C_k = 0.05;
% c1 = 0.056;  % courant coeff
% tmax = 1100; % total time
% C_k = 0.01;
% c1 = 0.0053;  % courant coeff
% tmax = 8000; % total time
dt = c1 * Re / N^2;
total_time = dt*tmax;
m = 2;
nut_inf = 0.31*(5-m)/(m+1)*sqrt(3-m)*C_k^(-3/2);
%% Initial condition
uk = zeros(1,N);
for k = 1:N   % uk(1) to uk(20)
    uk(k) = 1/k + 0i; %% default IC
%     uk(k) = 0 + 1/k*1i; %% IC2
end
%% solving uk for different modes
ut = zeros(1,N);
vt_star = zeros(1,N);
nu_t = zeros(1,N);
nu_eff = zeros(1,N);
err = zeros(1,tmax);
ut(1) = uk(1);    % ut(1) = uk(1) as u(1) is constant
tic
for t = 1:tmax
    for k = 2:N   % ut(2) to ut(20)
        up_iq_uq = 0i;
        for p = -N:N
            q = k - p;
            if q >= -N && q <= N
%             for q = -N:N
%                 if p + q == k
                    if p >= 1 && q >= 1
                        up_iq_uq = up_iq_uq + uk(p)*(0+1i)*q*uk(q);
                    elseif q < 0
                        up_iq_uq = up_iq_uq + uk(p)*(0+1i)*q*conj(uk(-q));
                    elseif p < 0
                        up_iq_uq = up_iq_uq + conj(uk(-p))*(0+1i)*q*uk(q);
                    end
%                 end
%             end
            end
        end
%         ut(k) = uk(k) + dt*( -up_iq_uq - k^2/Re*uk(k) ); % non-LES
        vt_star(k) = 1 + 34.5*exp(1)^(-3.03*(N/k));
        nu_t(k) = nut_inf*sqrt( uk(N)*conj(uk(N))/N )*vt_star(k);
        nu_eff(k) = 1/Re + nu_t(k); % for LES
%         nu_eff(k) = 1/Re; % for normal & DNS
        ut(k) = uk(k) + dt*( -up_iq_uq - k^2 * nu_eff(k) * uk(k) );
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
% semilogy(1:tmax,imag(sum_error))
xlim([0 tmax])
% ylim([1e-15 1e+0])
xlabel('time-step','interpreter','latex','fontsize',12)
ylabel('Error','interpreter','latex','fontsize',12)
grid on
legend_plot = sprintf('Euler, $N=20$ (LES), C$_k=0.4523$, $CFL$ = %.3f',c1);
% legend_plot = sprintf('Euler, $N=20$ (LES), C$_k=0.05$, $CFL$ = %.3f',c1);
% legend_plot = sprintf('Euler, $N=20$ (LES), C$_k=0.18$, $CFL$ = %.3f',c1);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Burgers equation: Error with increasing time-step','Interpreter','latex','fontsize',12)
print Burgers_ERROR-vs-Time-step.png -dpng
print Burgers_ERROR-vs-Time-step.eps -depsc
print Burgers_ERROR-vs-Time-step.pdf -dpdf
%%  Plot Error VS time graph
figure(2)
time_matrix = (dt:dt:total_time);
semilogy(time_matrix,real(sum_error))
% semilogy(time_matrix,imag(sum_error))
xlim([0 time_matrix(tmax)])
ylim([1e-15 1e+0])
xlabel('time','interpreter','latex','fontsize',12)
ylabel('Error','interpreter','latex','fontsize',12)
grid on
legend_plot = sprintf('Euler, $N=20$ (LES), C$_k=0.4523$,, $dt$ = %.4f',dt);
% legend_plot = sprintf('Euler, $N=20$ (LES), C$_k=0.05$, $dt$ = %.4f',dt);
% legend_plot = sprintf('Euler, $N=20$ (LES), C$_k=0.18$, $dt$ = %.4f',dt);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Burgers equation: Error with increasing time','Interpreter','latex','fontsize',12)
print Burgers_ERROR-vs-Time.png -dpng
print Burgers_ERROR-vs-Time.eps -depsc
print Burgers_ERROR-vs-Time.pdf -dpdf
%% Calculating & Plotting the real function 'u' over the real space
a_k = 2*real(ut); % Real parts from the Fourier coefficients
b_k = -2*imag(ut); % Imaginary parts from the Fourier coefficients
Nx = 7001; % Number of x values
u = zeros(1,Nx); % Initializing u
x_values = zeros(1,Nx); % x-values
for x = 1:Nx
u_old = 0;
    for k = 1:N
        x_values(x) = (x-1)/1000; % to start x-values from 0 to 150
        u(x) = u_old + a_k(k)*cos(k*(x-1)/1000) + b_k(k)*sin(k*(x-1)/1000);
        u_old = u(x);
    end
end
figure(4)
plot(x_values,u,'b-')
% plot(x_values,u,'-')
% plot(x_values,u,'m-')
xlim([0 2*pi])
set(gca,'xtick',(0:pi/2:2*pi)) % where to set the tick marks
set(gca,'xticklabels',{'0','\pi/2','\pi','3\pi/2','2\pi'}) % give them user-defined labels
ylim([-3 4])
legend_plot = sprintf('$u$, when $N=20$ (LES), $C_k=%.4f$',C_k);
legend(legend_plot,'Interpreter','latex','fontsize',11,'Location','Northeast')
title('Steady-state solution of Burgers equation','Interpreter','latex','fontsize',12)
xlabel('$x$','interpreter','latex','fontsize',12) 
ylabel('$u$','interpreter','latex','fontsize',12)
grid on
print Burgers_u.png -dpng
print Burgers_u.eps -depsc
print Burgers_u.pdf -dpdf
savefig('Burgers_u.fig')
hold on
%% solving Ek for different modes
E = zeros(1,N);
for k = 1:N
    E(k) = ut(k)*conj(ut(k));
end
%% visualize Ek for N = 20 (LES: Ck = 0.4523)
figure(3)
loglog(E,'b+-','LineWidth',1)
xlim([1 5*N])
ylim([1e-6 1e+0])
xlabel('$k$','interpreter','latex','fontsize',12) 
ylabel('$E_k$','interpreter','latex','fontsize',12)
grid on
legend_plot = sprintf('Euler, $N=20$, C$_k=0.4523$, $CFL$ = %.3f',c1);
% legend_plot = sprintf('Euler, $N=20$, C$_k=0.305$, $CFL$ = %.3f',c1);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Energy spectrum of Burgers equation''s steady-state solution','Interpreter','latex','fontsize',12)
hold on
%% visualize Ek for N = 20 (LES: Ck = 0.05)
% figure(3)
% % loglog(E,'g+-','LineWidth',1)
% loglog(E,'m+-','LineWidth',1)
% xlim([1 5*N])
% ylim([1e-6 1e+0])
% xlabel('$k$','interpreter','latex','fontsize',12) 
% ylabel('$E_k$','interpreter','latex','fontsize',12)
% grid on
% % legend_plot = sprintf('Euler, $N=20$, C$_k=0.05$, $CFL$ = %.3f',c1);
% legend_plot = sprintf('Euler, $N=20$, C$_k=0.01$, $CFL$ = %.3f',c1);
% legend(legend_plot,'Interpreter','latex','fontsize',11)
% title('Energy spectrum of Burgers equation''s steady-state solution','Interpreter','latex','fontsize',12)
% hold on
print Burgers.png -dpng
print Burgers.eps -depsc
print Burgers.pdf -dpdf
%% Calculating & Plotting the real function 'u' over the real space
% clear
% load ('data_Burgers_case-2')
% figure(3)
% loglog(E,'ok','LineWidth',1)
% xlim([1 1*N])
% ylim([1e-6 1e+0])
% % legend_plot = sprintf('Euler, $N=100$, $CFL$ = %.3f',c1);
% % legend(legend_plot,'Interpreter','latex','fontsize',11)
% xlabel('$k$','interpreter','latex','fontsize',12) 
% ylabel('$E_k$','interpreter','latex','fontsize',12)
% grid on
% print Burgers_E_comp.png -dpng
% print Burgers_E_comp.eps -depsc
% print Burgers_E_comp.pdf -dpdf
% savefig('Burgers_E_comp.fig')
% %% Plotting the Energy spectrum again to compare with DNS (N=100)
% figure(4)
% plot(x_values,u,'-k','LineWidth',1)
% xlim([0 2*pi])
% set(gca,'xtick',(0:pi/2:2*pi)) % where to set the tick marks
% set(gca,'xticklabels',{'0','\pi/2','\pi','3\pi/2','2\pi'}) % give them user-defined labels
% ylim([-3 4])
% % legend_plot = sprintf('$u$, when $N=100$',N);
% % legend(legend_plot,'Interpreter','latex','fontsize',11)
% xlabel('$x$','interpreter','latex','fontsize',12) 
% ylabel('$u$','interpreter','latex','fontsize',12)
% grid on
% print Burgers_u_comp.png -dpng
% print Burgers_u_comp.eps -depsc
% print Burgers_u_comp.pdf -dpdf
% savefig('Burgers_u.fig')