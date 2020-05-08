clear
close all

%load saved file from matcam.
out = load('hefwet + geometrie zonder excentriciteit.mat');
out2 = load('hefwet + geometrie.mat');
% 
% 
% roc_pitch = out.roc_pitch;
% roc_cam = out.roc_cam;
% deltaR = roc_pitch - roc_cam;
% 
% figure
% plot(out.theta, deltaR);
% xlabel('$\theta$ [rad]', 'Interpreter', 'LaTex');
% ylabel('$\rho_{pitch} - \rho_{cam}$ [$mm$]', 'Interpreter', 'LaTex');
% xlim([0,2*pi]);
% 
% max(out.pressure_angle);
% min(out.pressure_angle);
% verschil = abs(max(out.pressure_angle)-abs(min(out.pressure_angle)));

%% Bepalen van vermogen zonder excentriciteit
contactkracht = out.normalforce_tot;
drukhoek = out.pressure_angle;
hoeksnelheid = out.w;
xpitch = out.xpitch;
ypitch = out.ypitch;
nokstraal = sqrt(xpitch.^2+ypitch.^2).*10^(-3);

p = contactkracht.*sin(drukhoek).*nokstraal.*hoeksnelheid;
pp = ones(size(p));
pgem = pp.*mean(p);

figure('Name','Vermogen zonder excentriciteit')
plot(out.thetadegree,p,out.thetadegree,pgem)
title('Vermogen zonder excentriciteit');
xlabel('Nokhoek [°]');
ylabel('Vermogen [Watt]');
xlim([0,360]);
saveas(gcf,'vermogen.png');
ylim([-10^6,10^6]);

%% Bepalen van vermogen met excentriciteit
ck = out2.normalforce_tot;
dh = out2.pressure_angle;
hs = out2.w;
x2pitch = out2.xpitch;
y2pitch = out2.ypitch;
ns = sqrt(x2pitch.^2+y2pitch.^2).*10^(-3);

to1 = sqrt(ns.^2-(out2.exc*(10^-3))^2);
test = (tan(dh).*to1);
test2 = out2.exc*10^-3;
to2 = (tan(dh).*to1) - abs(out2.exc*10^(-3));

p2 = ck.*cos(dh).*hs.*to2;
pp2 = ones(size(p2));
p2gem = pp2.*mean(p2);

figure('Name','Vermogen met excentriciteit')
plot(out2.thetadegree,p2,out2.thetadegree,p2gem)
title('Vermogen met excentriciteit');
xlabel('Nokhoek [°]');
ylabel('Vermogen [Watt]');
xlim([0,360]);
saveas(gcf,'vermogenexc.png');
ylim([-10^6,10^6]);


%% Bepalen van verschil in vermogens
verschil_p = p-p2;
verschil_p_gem = pgem -p2gem;

figure('Name','Verschil in ogenblikkelijke vermogens')
plot(out2.thetadegree,verschil_p)
title('Verschil in ogenblikkelijke vermogens');
xlabel('Nokhoek [°]');
ylabel('Vermogen [Watt]');
xlim([0,360]);
saveas(gcf,'delta_vermogen.png');

figure('Name','Verschil in gemiddelde vermogens')
plot(out2.thetadegree,verschil_p_gem)
title('Verschil in gemiddelde vermogens');
xlabel('Nokhoek [°]');
ylabel('Vermogen [Watt]');
xlim([0,360]);
saveas(gcf,'delta_gem_vermogen.png');
ylim([-10^-13,10^-13]);


% %% Bepalen vliegwiel
% koppel = (p)./hoeksnelheid;
% kgem = (pgem./hoeksnelheid).*ones(1,length(out.theta));
% 
% figure
% plot(out.theta,koppel,out.theta,kgem)
% title('Koppel en het gemiddeld koppel');
% xlabel('Nokhoek [rad]');
% ylabel('Koppel [Nm]');
% saveas(gcf,'koppel.png');
% 
% wnom=out.w;
% w_min = 0.95*wnom;
% w_max = 1.05*wnom;
% K = (w_max-w_min)/wnom;
% 
% A = cumtrapz(out.theta,koppel-kgem);
% figure
% plot(out.theta,A)
% xlabel('Nokhoek [rad]');
% ylabel('Koppel [Nm]');
% 
% [~,theta_max] = max(A);
% theta_MD = out.thetadegree(theta_max);
% [~,theta_min] = min(A);
% theta_mD = out.thetadegree(theta_min);
% 
% [km,KM] = min(koppel);
% 
% figure
% plot(out.theta,koppel,out.theta,kgem,out.theta(theta_max),koppel(theta_max),'*r',out.theta(theta_min),koppel(theta_min),'*r',out.theta(KM),koppel(KM),'*r')
% title('Koppel en het gemiddeld koppel');
% xlabel('Nokhoek [rad]');
% ylabel('Koppel [Nm]');
% saveas(gcf,'koppelben.png');
% 
% integrand = koppel-kgem;
% Amax = trapz(out.theta(theta_min:theta_max),integrand(theta_min:theta_max));
% I_vliegwiel = Amax/(K*wnom^2);
% 
% oppervlakte = (out.theta(theta_max)-out.theta(theta_min))*abs((km-kgem(1)))/2;
% I_controle = oppervlakte/(K*wnom^2);
% straal = ((2*I_vliegwiel)/(pi*0.05*7800))^(1/4);
% 
% 
% 
% %% Deel3: numerieke simulatie single-rise
%  zeta = 0.08;
%  lambda = 9.38;
%  t_begin=7000;
%  t_end=18000;
% 
% teller = (2*pi*lambda)^2;
% noemer = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
% 
% tau = (out2.theta(14001:36000)-out2.theta(14001))/(70*2*pi/360);
% theta = out2.S(14001:36000)/25;
% theta0 = 1;
% theta_dot0 = 0;
% [A,B,C,D] = tf2ss(teller,noemer);
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