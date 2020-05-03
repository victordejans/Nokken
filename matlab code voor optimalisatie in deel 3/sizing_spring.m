function k = sizing_spring
%returns the ideal spring constant k when given a pre load Fv0

close all
clear

%defining variables and parameters
theta=(1:1:360)';
S = hefwet;
Ffunc = external_forces;
Fv0 = 00*ones(size(theta));
m = 20;
omega = 1;
acc = acceleration(S);


%calculating the function that needs to be maximized
F = (-Ffunc-Fv0-m*omega^2*acc)./S;



k = max(F);


plot(theta,F)


end