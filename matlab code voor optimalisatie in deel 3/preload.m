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