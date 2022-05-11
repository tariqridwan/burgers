%% Burgers version 5.O: steady-state calculation_Optimized
close all
clear
%% Physical description
Re = 40; % Reynolds number
% N = 20; % total number of fourier modes
% c1 = 0.3;  % courant coeff
% tmax = 900; % total time
N = 100; % total number of fourier modes
c1 = 0.76;  % courant coeff
tmax = 1750; % total time
dt = c1 * Re / N^2;
total_time = dt*tmax;
%% Initial condition
uk = zeros(1,N);
for k = 1:N   % uk(1) to uk(20)
    uk(k) = 1/k + 0i; %% default IC
%     uk(k) = 0 + 1/k*1i; %% IC2
end
%% solving uk for different modes
ut = zeros(1,N);
% err = zeros(1,N);
% err = zeros(N,tmax);
err = zeros(1,tmax);
ut(1) = uk(1);      % ut(1) = uk(1) as u(1) is constant
tic
for t = 1:tmax
    for k = 2:N  % ut(2) to ut(N)
        up_iq_uq = 0i;  %Initializing the term
        for p = -N:N  % setting q values: -N < p < N
            q = k - p;  % Instead of q = -N:N
            if q >= -N && q <= N  % limiting q values: -N < q < N
                    if p >= 1 && q >= 1  % when p & q are positive
                        up_iq_uq = up_iq_uq + uk(p)*(0+1i)*q*uk(q);
                    elseif q < 0  % when only q is negative (p is positive)
                        up_iq_uq = up_iq_uq + uk(p)*(0+1i)*q*conj(uk(-q));
                    elseif p < 0  % when only p is negative (q is positive)
                        up_iq_uq = up_iq_uq + conj(uk(-p))*(0+1i)*q*uk(q);
                    end
            end
        end        
        ut(k) = uk(k) + dt*( -up_iq_uq - k^2/Re*uk(k) );  % calculate u
        err(t) = err(t) + ( ut(k) - uk(k) )^2;  % calculate conv. error
        uk(k) = ut(k);  % set u for next time-step
    end
end
toc
CPUtime = toc
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
% legend({'Error with increasing time-step'},'Interpreter','latex','fontsize',10,'Location','Northeast')
% legend_plot = sprintf('Euler, $N=20$, $CFL$ = %.3f',c1);
legend_plot = sprintf('Euler, $N=100$, $CFL$ = %.3f',c1);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Burgers equation: Error with increasing time-step','Interpreter','latex','fontsize',12)
% title_plot = sprintf('$N$ = %d',N);
% title(title_plot,'Interpreter','latex','fontsize',12)
print Burgers_ERROR-vs-Time-step.png -dpng
print Burgers_ERROR-vs-Time-step.eps -depsc
print Burgers_ERROR-vs-Time-step.pdf -dpdf
savefig('Burgers_ERROR-vs-Time-step.fig')
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
% legend({'Error with increasing time'},'Interpreter','latex','fontsize',10,'Location','Northeast')
% legend_plot = sprintf('Euler, $N=20$, $dt$ = %.4f',dt);
legend_plot = sprintf('Euler, $N=100$, $dt$ = %.4f',dt);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Burgers equation: Error with increasing time','Interpreter','latex','fontsize',12)
% title_plot = sprintf('$N$ = %d',N);
% title(title_plot,'Interpreter','latex','fontsize',12)
print Burgers_ERROR-vs-Time.png -dpng
print Burgers_ERROR-vs-Time.eps -depsc
print Burgers_ERROR-vs-Time.pdf -dpdf
savefig('Burgers_ERROR-vs-Time.fig')
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
% legend_plot = sprintf('Euler, $N=100$, $CFL$ = %.3f',c1);
% legend(legend_plot,'Interpreter','latex','fontsize',11)
% title('Energy spectrum of Burgers equation''s steady-state solution','Interpreter','latex','fontsize',12)
% % title_plot = sprintf('$CFL$ = %.3f',c1);
% % title(title_plot,'Interpreter','latex','fontsize',12)
% print Burgers.png -dpng
% print Burgers.eps -depsc
% print Burgers.pdf -dpdf
% savefig('Burgers.fig')
%% visualize Ek for N = 20
figure(3)
loglog(E,'+r-','LineWidth',1)
xlim([1 5*N])
ylim([1e-6 1e+0])
xlabel('$k$','interpreter','latex','fontsize',12) 
ylabel('$E_k$','interpreter','latex','fontsize',12)
grid on
% legend({'$N=20$'},'Interpreter','latex','fontsize',10,'Location','Northeast')
legend_plot = sprintf('Euler, $N=20$, $CFL$ = %.3f',c1);
legend(legend_plot,'Interpreter','latex','fontsize',11)
title('Energy spectrum of Burgers equation''s steady-state solution','Interpreter','latex','fontsize',12)
% title_plot = sprintf('$CFL$ = %.3f',c1);
% title(title_plot,'Interpreter','latex','fontsize',12)
print Burgers.png -dpng
print Burgers.eps -depsc
print Burgers.pdf -dpdf
savefig('Burgers.fig')
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
plot(x_values,u,'-')
xlim([0 2*pi])
set(gca,'xtick',(0:pi/2:2*pi)) % where to set the tick marks
set(gca,'xticklabels',{'0','\pi/2','\pi','3\pi/2','2\pi'}) % give them user-defined labels
ylim([-3 4])
legend_plot = sprintf('$u$, when $N=%d$',N);
legend(legend_plot,'Interpreter','latex','fontsize',11,'Location','Northwest')
title('Steady-state solution of Burgers equation','Interpreter','latex','fontsize',12)
xlabel('$x$','interpreter','latex','fontsize',12) 
ylabel('$u$','interpreter','latex','fontsize',12)
grid on
print Burgers_u.png -dpng
print Burgers_u.eps -depsc
print Burgers_u.pdf -dpdf
savefig('Burgers_u.fig')