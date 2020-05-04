function optimal_spring

%reading data from .mat file
out=load('hefwet + geometrie.mat');
alpha_=out.pressure_angle;
for i = 1:length(alpha_)/100
    alpha(i)=alpha_(100*i)*180/pi;
end
omega=out.w;
S_=out.S;
for i = 1:length(S_)/100
    S(i)=S_(100*i);
end
acc_=out.A;
for i = 1:length(acc_)/100
    acc(i)=acc_(100*i);
end
Ffunc_=out.extload;
for i = 1:length(Ffunc_)/100
    Ffunc(i)=Ffunc_(100*i);
end
m=20;
syms k; %de te zoeken veerconstante
syms 'Fv0'; %de te zoeken preload

% 
% Ffunc(1)
% S(1)
% acc(1)
% cos(alpha(1))
% (Ffunc(1)+S(1)*k+Fv0+m*acc(1))/cos(alpha(1))
N = @(Fv0,k)(Ffunc+Fv0+k*S+m*omega^2*acc)./cos(alpha);

% hallo =
% figure
% plot(hallo)
integral = @(Fv0,k)sum(N)/360;

    


