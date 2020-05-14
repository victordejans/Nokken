clear
close all

%load saved file from matcam.
out = load('hefwet + geometrie zonder excentriciteit.mat');
out2 = load('hefwet + geometrie.mat');


%% Determining k_f
 m = out.mass; % mass of 20 kg
 zeta = 0.054;
 omega = out.w;
 
 lambda_min = 0.75/zeta;  %lambda has to be greater than 13.8889
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

T = 36000;  % 1 cycle consists of 36000 data points

t_begin = T*theta_end/(2*pi);  % start of the free response
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

theta_equi_begin = 0;  % in degrees for ease of use
theta_equi_end = 150;  % also in degrees
cycles = 15;  % amount of cycles after which we certainly have achieved equilibrium
t_multi_begin = theta_equi_begin*100 + 36000*cycles;
t_multi_end = theta_equi_end*100 + 36000*cycles;

numerator_multi = (2*pi*lambda_multi)^2;
denominator_multi = [1, 2*zeta*(2*pi*lambda_multi), (2*pi*lambda_multi)^2];
sys_multi = tf(numerator_multi, denominator_multi);

tau_multi = linspace(0,25,36000*25);  % extended tau so we can also study steady-state
theta_multi= repmat(out2.S,[1,25])./35;  % extended normalized lift for steady-state analysis
gamma_multi= lsim(sys_multi,theta_multi,tau_multi)';

figure('Name','Response gamma Multirise','Position', [100, 100, 500, 500])

subplot(2,1,1)
plot(tau_multi(t_multi_begin:t_multi_end),theta_multi(t_multi_begin:t_multi_end),'LineWidth',1);
hold on
plot(tau_multi(t_multi_begin:t_multi_end),gamma_multi(t_multi_begin:t_multi_end),'LineWidth',1);
title('Excitation and response')
xlabel('cycle')
ylabel('Dimensionless lift [-]')
axis('tight');
legend('\theta','\gamma')

subplot(2,1,2)
plot(tau_multi(t_multi_begin:t_multi_end),gamma_multi(t_multi_begin:t_multi_end) - theta_multi(t_multi_begin:t_multi_end),'LineWidth',1);
title('Difference of excitation and response')
xlabel('cycle')
ylabel('\gamma - \theta  [-]')
axis('tight');
saveas(gcf,'multirise.png');


%% Comparing multirise to single rise

% gamma_norm = gamma(1:multi_end-multi_begin+1);

figure('Name','Comparing Single and Multi rise','Position', [100, 100, 500, 500])

subplot(2,1,1)
plot(tau_multi(t_multi_begin:t_multi_end),gamma_multi(t_multi_begin:t_multi_end),'LineWidth',1);
hold on
plot(tau_multi(t_multi_begin + 1:t_multi_end),gamma(1:t_end),'LineWidth',1);
title('Multi and Single rise')
xlabel('\tau [-]')
ylabel('Dimensionless lift [-]')
axis ('tight');
legend('multi rise','single rise')


subplot(2,1,2)
plot(tau_multi(t_multi_begin + 1:t_multi_end), gamma_multi(t_multi_begin + 1:t_multi_end)' - gamma(1:t_end),'LineWidth',1)
title('Difference of Multi and Single rise')
xlabel('\tau [-]')
ylabel('\gamma_{multi} - \gamma_{single} [-]')
axis ('tight');
saveas(gcf,'multirisesinglerisevergelijken.png');

theta_deg_free_start = 100; % Angle in degrees for when free response approximately starts
t_multi_free_begin = theta_deg_free_start*100 + 36000*cycles;

figure('Name','Comparing Single and Multi rise free response','Position', [100, 100, 500, 500])

subplot(2,1,1)
plot(tau_multi(t_multi_free_begin:t_multi_end),gamma_multi(t_multi_free_begin:t_multi_end),'LineWidth',1);
hold on
plot(tau_multi(t_multi_free_begin:t_multi_end),gamma(t_begin:t_end),'LineWidth',1);
title('Multi and Single rise free response')
xlabel('\tau [-]')
ylabel('Dimensionless lift [-]')
axis ('tight');
legend('multi-rise','single-rise')


subplot(2,1,2)
plot(tau_multi(t_multi_free_begin:t_multi_end), gamma_multi(t_multi_free_begin:t_multi_end)' - gamma(t_begin:t_end),'LineWidth',1)
title('Difference in Mr and Sr free response')
xlabel('\tau [-]')
ylabel('\gamma_{m} - \gamma_{s} [-]')
axis ('tight');
saveas(gcf,'multirisesinglerisevrijeresponsvergelijken.png');
 
%% Final contact force
gamma_ref = gamma_multi(t_multi_begin + 1:t_multi_begin + 36000);  % reference gamma for our other properties
theta_contact = 2*pi./length(gamma_ref); % artificial theta
S_contact = gamma_ref.*35; %Accurate, actual displacement
V_contact = gradient(S_contact,theta_contact); % actual speed
A_contact = gradient(V_contact,theta_contact); % actual acceleration

N_load = out2.extload./cos(out2.pressure_angle);   
N_spring = (S_contact*out2.springconstant + out2.springpreload)./cos(out2.pressure_angle);             
N_acc = (m*(A_contact./(10^6))*(omega^2))./cos(out2.pressure_angle);
N = N_load + N_spring + N_acc;

figure('Name', 'Final contact force')

subplot(2,1,1)
plot(out2.theta, N,'LineWidth',1)
hold on
plot(out2.theta, out2.normalforce_tot,'LineWidth',1)
hold on
xlabel('\theta [rad]')
ylabel('Contact force N [N]')
axis ([0 2*pi 0 800]);
legend('N for deformable body','N for rigid body')

subplot(2,1,2)
plot(out2.theta, N - out2.normalforce_tot,'LineWidth',1)
xlabel('\theta [rad]')
ylabel('Difference in contact force N [N]')
axis ('tight');

saveas(gcf,'Contactkracht.png');
