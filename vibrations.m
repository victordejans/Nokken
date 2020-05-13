clear
close all

%load saved file from matcam.
out = load('hefwet + geometrie zonder excentriciteit.mat');
out2 = load('hefwet + geometrie.mat');


%% Determining k_f
 m = out.mass; % mass of 20 kg
 zeta = 0.054;
 omega = out.w;
 
 lambda_min = 0.75/zeta  %lambda has to be greater than 13.8889
 lambda = 14;  % So we choose lambda to be 14
 
 % define the critical lift
 theta_start = 0;
 theta_end = 100 * pi/180;  %angles in radians
 S_start = 0;
 S_end = 35; % lift in mm
 
 t_1 = (theta_end - theta_start)/omega;
 t_n = t_1/lambda;
 omega_n = 2*pi/t_n;
 k_f = ((omega_n)^2)*m  %result is a spring constant of about 2*10^6 N/mm
 
 
%% Numerical single rise analysis

T = 36000;  % 1 cycle consists of 36000 data points

t_begin = T*theta_end/(2*pi)  % start of the 'vrije respons'
t_end = T*150/360;  % start of the next rise

tau = out2.theta(1:t_end)/theta_end;
theta = out2.S(1:t_end)/S_end;

numerator = (2*pi*lambda)^2;
denominator = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(numerator, denominator);
gamma = lsim(sys,theta,tau);


figure('Name','Response gamma', 'Position', [100, 100, 500, 800])

subplot(3,1,1);
plot(tau,theta,'LineWidth',1);
hold on
plot(tau,gamma,'LineWidth',1);
title('Excitation and response as a function of \tau')
xlabel('\tau [-]')
ylabel('Dimensionless rise [-]')
axis('tight');
legend('\theta','\gamma')

subplot(3,1,2);
plot(tau,gamma'-theta,'LineWidth',1);
axis('tight');
title('Difference of response and excitation')
xlim([0 tau(end)]);
xlabel('\tau [-]')
ylabel('\gamma - \theta [-]')

subplot(3,1,3)
plot(tau(t_begin:t_end),theta(t_begin:t_end),'LineWidth',1);
hold on
plot(tau(t_begin:t_end),gamma(t_begin:t_end),'LineWidth',1);
title('Free response: excitation and response')
xlabel('\tau [-]')
ylabel('Dimensionless rise [-]')
axis('tight');
legend('\theta','\gamma')
xlim([tau(t_begin) tau(t_end)])

saveas(gcf,'SingleriseNumeriek.png');


%% Approximate single rise analysis

N=3;
Q=(2*pi)^2;
A_1 = Q/(2*pi*lambda).^N;  %approximate A_1
bounds_exp = A_1*exp(-zeta*2*pi*lambda*(tau-1));

figure('Name','Approximate analysis', 'Position', [100, 100, 500, 800])

subplot(2,1,1)
plot(tau(t_begin:t_end),bounds_exp(t_begin:t_end),'LineWidth',1);
hold on
plot(tau(t_begin:t_end),gamma(t_begin:t_end) - 1,'LineWidth',1);
xlim([1 1.5])
plot(tau(t_begin:t_end),theta(t_begin:t_end) - 1,'LineWidth',1)
title('bounding exponential')
xlabel('\tau [-]')
ylabel('Dimensionless rise [-]')
legend('exp_{bounds}','\gamma','\theta')

subplot(2,1,2)
plot(tau(t_begin:t_end), bounds_exp(t_begin:t_end) + 1 - gamma(t_begin:t_end)','LineWidth',1);
hold on
xlim([1 1.5])
title('Approximation error')
xlabel('\tau [-]')
ylabel('exp_{bounds} - \gamma [-]')
axis('tight');
saveas(gcf,'SingleriseBenadering.png');


%% Multi rise analysis
lambda_multi = 1/t_n;  %We need a new lambda for multi rise analysis

numerator_multi = (2*pi*lambda_multi)^2;
denominator_multi = [1, 2*zeta*(2*pi*lambda_multi), (2*pi*lambda_multi)^2];
sys_multi = tf(numerator_multi, denominator_multi);

tau_multi = out2.theta/(2*pi);  % normalized time over 1 whole cycle
theta_multi= out2.S/35;  % normalzsed lift over 1 whole cycle
gamma_multi= lsim(sys_multi,theta_multi,tau_multi)';

figure('Name','Response gamma Multirise','Position', [100, 100, 500, 500])

subplot(2,1,1)
plot(tau_multi,theta_multi,'LineWidth',1);
hold on
plot(tau_multi,gamma_multi,'LineWidth',1);
title('Excitation and response')
xlabel('cycle')
ylabel('Dimensionless lift [-]')
legend('\theta','\gamma')

subplot(2,1,2)
plot(tau_multi,gamma_multi - theta_multi,'LineWidth',1);
title('Difference of excitation and response')
xlabel('cycle')
ylabel('\gamma - \theta  [-]')
saveas(gcf,'multirise.png');

% 
% %% vergelijken multirise - singlerise
%  begin_angle = 140;
%  end_angle = 250;
%  multi_begin = 360*100*20+begin_angle*100;
%  multi_end = 360*100*20+end_angle*100;
%  gamma_norm = gamma(1:multi_end-multi_begin+1);
% 
% figure('Name','Vergelijking Single- en Multi-rise','Position', [100, 100, 500, 500])
% subplot(2,1,1)
% plot(tau_mr(multi_begin:multi_end),gamma_mr(multi_begin:multi_end),'LineWidth',1);
% hold on
% plot(tau_mr(multi_begin:multi_end),gamma_norm,'LineWidth',1);
% title('Multi- en Single-rise')
% xlabel('\tau [-]')
% ylabel('Heffing [-]')
% legend('multi-rise','single-rise')
% axis tight;
% subplot(2,1,2)
% plot(tau_mr(multi_begin:multi_end), gamma_mr(multi_begin:multi_end)' - gamma_norm,'LineWidth',1)
% title('Verschil Multi- en Single-rise')
% xlabel('\tau [-]')
% ylabel('\gamma_{multi} - \gamma_{single} [-]')
% axis tight;
% saveas(gcf,'vglmultirisesinglerise.png');
% 
%  begin_angle_new = 210;
%  multi_begin_new = 360*100*20+begin_angle_new*100;
%  multi_delta = multi_end - multi_begin;
%  gamma_norm = gamma((begin_angle_new-begin_angle)/(end_angle-begin_angle)*multi_delta:multi_delta);
% 
% figure('Name','Vergelijking Single- en Multi-rise','Position', [100, 100, 500, 500])
% subplot(2,1,1)
% plot(tau_mr(multi_begin_new:multi_end),gamma_mr(multi_begin_new:multi_end),'LineWidth',1);
% hold on
% plot(tau_mr(multi_begin_new:multi_end),gamma_norm,'LineWidth',1);
% title('Multi- en Single-rise: vrije respons')
% xlabel('\tau [-]')
% ylabel('Heffing [-]')
% legend('multi-rise','single-rise')
% axis tight;
% subplot(2,1,2)
% plot(tau_mr(multi_begin_new:multi_end), gamma_mr(multi_begin_new:multi_end)' - gamma_norm,'LineWidth',1)
% title('Verschil Multi- en Single-rise: vrije respons')
% xlabel('\tau [-]')
% ylabel('\gamma_{m} - \gamma_{s} [-]')
% axis tight;
% saveas(gcf,'vglmultirisesinglerisevrijerespons.png');
%  
% %% nieuwe contactkracht
% gamma_analyse=gamma_mr(855008:891007);
% h=2*pi./length(gamma_analyse); % rad
% S_mr=gamma_analyse.*40; %staat in mm
% V_mr=gradient(S_mr,h); %mm/rad
% A_mr=gradient(V_mr,h); %mm/rad^2
% normalforce_spring = (out2.S*out2.springconstant+out2.springpreload)./cos(out2.pressure_angle);
% normalforce_load = out2.extload./cos(out2.pressure_angle);                
% normalforce_acc = out2.mass*A_mr./1000*(out2.w^2)./cos(out2.pressure_angle);
% normalforce_tot = normalforce_spring + normalforce_load + normalforce_acc;
% figure('Name', 'oud vs. new')
% plot(out2.theta,normalforce_tot,'LineWidth',1)
% hold on
% plot(out2.theta,out2.normalforce_tot,'LineWidth',1)
% %kleiner dan nul --> groter veertje
% hold on
% plot(out2.theta,zeros(size(out2.theta)))
% xlabel('\theta [rad]')
% ylabel('Normaalcontactkracht N [N]')
% legend('Vervormbare volger','Onvervormbare volger')
% xlim([0 2*pi])
% ylim([-100 500])
% saveas(gcf,'normaalkrachten.png');
% 
%  