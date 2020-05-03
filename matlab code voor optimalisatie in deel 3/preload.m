function Fv0 = preload

close all
clear

%defining variables and parameters
theta=(1:1:360)';
S = hefwet;
Ffunc = external_forces;
Fv0 = 0*ones(size(theta));
m = 20;
omega = 1;
acc = acceleration(S);
alpha = pressure_angle(S);
k=28.4;

Fv0_min = (-Ffunc-m*omega^2*acc-k*S)./cos(alpha*pi/180);

figure
plot(theta,Fv0_min)

Fv0 = max(Fv0_min);

end