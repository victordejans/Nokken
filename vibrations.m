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
 k_f = ((omega_n)^2)*m;  %result is a spring constant of about 2*10^6 N/mm
 
 
%% Numerical single rise analysis

%  t_begin=7000;
%  t_end=18000;

numerator = (2*pi*lambda)^2;
denominator = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];

% tau = (out2.theta(14001:36000)-out2.theta(14001))/(70*2*pi/360);
% theta = out2.S(14001:36000)/25;
% theta0 = 1;
% theta_dot0 = 0;
 [A,B,C,D] = tf2ss(numerator,denominator);
% X0 = [1/C(2)*theta_dot0; 1/C(2)*theta0];
% lsim(A,B,C,D, theta, tau, X0)
% gamma = lsim(A,B,C,D, theta, tau, X0);
% 
% figure('Name','Respons gamma', 'Position', [100, 100, 500, 800])
% subplot(3,1,1);
% plot(tau,theta,'LineWidth',1);
% hold on
% plot(tau,gamma,'LineWidth',1);
% title('Excitatie en respons ifv \tau')
% xlabel('\tau [-]')
% ylabel('Heffing [-]')
% axis('tight');
% legend('\theta','\gamma')
% subplot(3,1,2);
% plot(tau,gamma'-theta,'LineWidth',1);
% axis('tight');
% title('Verschil respons en excitatie ifv \tau')
% xlim([0 tau(end)]);
% xlabel('\tau [-]')
% ylabel('\gamma-\theta [-]')
% subplot(3,1,3)
% plot(tau(t_begin:t_end),theta(t_begin:t_end),'LineWidth',1);
% hold on
% plot(tau(t_begin:t_end),gamma(t_begin:t_end),'LineWidth',1);
% title('Vrije respons: excitatie en respons ifv \tau')
% xlabel('\tau [-]')
% ylabel('Heffing [-]')
% axis('tight');
% legend('\theta','\gamma')
% xlim([tau(t_begin) tau(t_end)])
% saveas(gcf,'SingleriseNumeriek.png');
% 
% 
% %% Benaderende analyse single-rise
% N=3;
% Q=(2*pi)^2;
% A1_ben=Q/(2*pi*lambda).^N;
% exp_plot=A1_ben.*exp(-zeta*2*pi.*lambda.*(tau-1));
% 
% figure('Name','Benaderend', 'Position', [100, 100, 500, 800])
% subplot(2,1,1)
% plot(tau(t_begin:t_end),exp_plot(t_begin:t_end),'LineWidth',1);
% hold on
% plot(tau(t_begin:t_end),gamma(t_begin:t_end),'LineWidth',1);
% xlim([1 2])
% plot(tau,zeros(size(tau)),'LineWidth',1)
% title('Exponentiële omhullende')
% xlabel('\tau [-]')
% ylabel('Heffing [-]')
% legend('exp_{benaderend}','\gamma','\theta')
% subplot(2,1,2)
% plot(tau(t_begin:t_end), exp_plot(t_begin:t_end) - gamma(t_begin:t_end)','LineWidth',1);
% hold on
% xlim([1 2])
% title('Verschil benadering en exacte methode')
% xlabel('\tau [-]')
% ylabel('exp_{benaderend} - \gamma [-]')
% axis('tight');
% saveas(gcf,'SingleriseBenadering.png');
% 
% 
% %% Multi-rise analyse
%  lambda_accent = 48.23;
% 
% teller = (2*pi*lambda_accent)^2;
% noemer = [1, 2*zeta*(2*pi*lambda_accent), (2*pi*lambda_accent)^2];
% sys = tf(teller, noemer);
% 
% tau_mr = linspace(0,25,36000*25);
% theta_mr=repmat(out.S,[1,25])./25;
% gamma_mr=lsim(sys,theta_mr,tau_mr)';
% 
% figure('Name','Respons gamma Multirise','Position', [100, 100, 500, 500])
% subplot(2,1,1)
% plot(tau_mr,theta_mr,'LineWidth',1);
% hold on
% plot(tau_mr,gamma_mr,'LineWidth',1);
% axis([20 21 -0.1 1.1])
% title('Excitatie en respons')
% xlabel('periode')
% ylabel('Heffing [-]')
% legend('\theta','\gamma')
% subplot(2,1,2)
% plot(tau_mr,gamma_mr-theta_mr,'LineWidth',1);
% title('Verschil excitatie en respons')
% xlabel('periode')
% ylabel('\gamma-\theta  [-]')
% axis([20 21 -0.007 0.007])
% saveas(gcf,'multirise.png');
% 
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