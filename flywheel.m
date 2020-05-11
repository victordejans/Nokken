function flywheel

close all; clear;

% Load our cam data
cam = load('hefwet + geometrie.mat');
N = cam.normalforce_tot;
omega = cam.w;
S = cam.S*0.001;
radians = cam.theta;
degrees = cam.theta/pi*180;
P = power_cam;
t = size(degrees,2);

%Instantaneous torque
M_inst = P ./ omega;

%Average torque
M_av = sum(M_inst)/t;

figure
plot(degrees,M_inst)
hold on
yline(M_av)
title("Instantaneous and average torque (Nm)")
xlabel("Cam angle (degrees)")
ylabel("Torque (Nm)")

%Work surplus
A = cumtrapz(radians,M_inst-M_av);
figure
plot(degrees,A)
xlabel('Cam angle (degrees)');
ylabel('Work surplus (J)');

%Theta_min and theta_max
[~,theta_min] = min(A); %index in theta array with minimal A
[~,theta_max] = max(A); %index in theta array with maximal A

%Maximal work surplus
A_max = trapz(radians(theta_min:theta_max),M_inst(theta_min:theta_max)-M_av);

%Fluctuation coefficient
omega_min = 0.95*omega;
omega_max = 1.05*omega;
K = (omega_max-omega_min)/omega;

%Moment of inertia flywheel
I = A_max/(K*omega^2);
I





end