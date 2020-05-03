function alpha = pressure_angle(S)

close all

t = size(S);
theta = (1:1:t)';
vel = velocity(S);
R0 = 60;
e = 4.5;

vel(2)


alpha = atan((vel-e)./(sqrt(R0^2-e^2)+S))*180/pi;


figure
plot(theta,alpha)


end
