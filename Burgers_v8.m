% This program solves the Burgers equation using Fourier analysis
% and explicit Runge-Kutta 2nd order method
% with Large Eddy Simulation (LES) technique
%            -------
% written by Tariq Ridwan: https://tariqridwan.github.io/
% Barcelona Supercomputing Center // Universitat Politècnica de Catalunya
%% Burgers version 7: RK2 steady-state calculation
close all
clear
%% Physical description
Re = 40; % Reynolds number
N = 20; % total number of fourier modes
C_k = 0.4523;
c1 = 0.201;  % courant coeff
tmax = 430; % total time
% C_k = 0.05;
% c1 = 0.116;  % courant coeff
% tmax = 550; % total time
dt = c1 * Re / N^2;
total_time = dt*tmax;
m = 2;
nut_inf = 0.31*(5-m)/(m+1)*sqrt(3-m)*C_k^(-3/2);
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
vt_star = zeros(1,N);
nu_t = zeros(1,N);
nu_t_0 = zeros(1,N);
nu_t_1 = zeros(1,N);
nu_eff = zeros(1,N);
nu_eff_0 = zeros(1,N);
nu_eff_1 = zeros(1,N);
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
        vt_star(k) = 1 + 34.5*exp(1)^(-3.03*(N/k));
        nu_t_0(k) = nut_inf*sqrt( uk(N)*conj(uk(N))/N )*vt_star(k);
        nu_eff_0(k) = 1/Re + nu_t_0(k); % for LES
%         nu_eff_0(k) = 1/Re; % for normal & DNS
        Ri_0(k) = -up_iq_uq - k^2 * nu_eff_0(k) * uk(k);
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
%         vt_star(k) = 1 + 34.5*exp(1)^(-3.03*(N/k)); % NO NEED TO REPEAT
        nu_t_1(k) = nut_inf*sqrt( ui_1(N)*conj(ui_1(N))/N )*vt_star(k);
        nu_eff_1(k) = 1/Re + nu_t_1(k); % for LES
%         nu_eff_1(k) = 1/Re; % for normal & DNS
        Ri_1(k) = -up_iq_uq - k^2 * nu_eff_1(k) * ui_1(k); % RK 2nd stage ends
        
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
% ylim([1e-15 1e+0])
xlabel('time-step','interpreter','latex','fontsize',12)
ylabel('Error','interpreter','latex','fontsize',12)
grid on
% legend({'Error with increasing time-step'},'Interpreter','latex','fontsize',10,'Location','Northeast')
legend_plot = sprintf('RK2, $N=20$ (LES), C$_k=0.4523$, $CFL$ = %.3f',c1);
% legend_plot = sprintf('RK2, $N=20$ (LES), C$_k=0.05$, $CFL$ = %.3f',c1);
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
legend_plot = sprintf('RK2, $N=20$ (LES), C$_k=0.4523$, $dt$ = %.4f',dt);
% legend_plot = sprintf('RK2, $N=20$ (LES), C$_k=0.05$, $dt$ = %.4f',dt);
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
%% visualize Ek for N = 20 (LES: Ck = 0.4523)
figure(3)
loglog(E,'b+-','LineWidth',1)
xlim([1 5*N])
ylim([1e-6 1e+0])
xlabel('$k$','interpreter','latex','fontsize',12) 
ylabel('$E_k$','interpreter','latex','fontsize',12)
grid on
% legend({'$N=20$'},'Interpreter','latex','fontsize',10,'Location','Northeast')
legend_plot = sprintf('RK2, $N=20$ (LES), C$_k=0.4523$, $CFL$ = %.3f',c1);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Energy spectrum of Burgers equation''s steady-state solution','Interpreter','latex','fontsize',12)
% title_plot = sprintf('$CFL$ = %.3f',c1);
% title(title_plot,'Interpreter','latex','fontsize',12)
print Burgers.png -dpng
print Burgers.eps -depsc
print Burgers.pdf -dpdf
%% visualize Ek for N = 20 (LES: Ck = 0.05)
% figure(3)
% loglog(E,'g+-','LineWidth',1)
% xlim([1 5*N])
% ylim([1e-6 1e+0])
% xlabel('$k$','interpreter','latex','fontsize',12) 
% ylabel('$E_k$','interpreter','latex','fontsize',12)
% grid on
% % legend({'$N=20$'},'Interpreter','latex','fontsize',10,'Location','Northeast')
% legend_plot = sprintf('Euler, $N=20$, C$_k=0.05$, $CFL$ = %.3f',c1);
% legend(legend_plot,'Interpreter','latex','fontsize',11)
% title('Energy spectrum of Burgers equation''s steady-state solution','Interpreter','latex','fontsize',12)
% % title_plot = sprintf('$CFL$ = %.3f',c1);
% % title(title_plot,'Interpreter','latex','fontsize',12)
% print Burgers.png -dpng
% print Burgers.eps -depsc
% print Burgers.pdf -dpdf