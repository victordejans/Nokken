function alpha = pressure_angle(S)

t = size(S);
theta = (1:1:t)';
vel = velocity(S);
R0 = 60;
e = 4.5;

alpha = zeros(t);

for i = 1:t
    alpha(i) = atan((vel(i)-e)/(sqrt(R0^2-e^2)+S(i)))*180/pi;
    i
%     if i == 13
%         vel(i)
%         e
%         R0
%         S(i)
%     end
end



figure
plot(theta,alpha)


end
