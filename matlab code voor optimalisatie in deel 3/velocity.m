function vel = velocity(S)

t = size(S);
theta = (1:1:t);
S = [S;S(1)];

vel = zeros(t);
for i = 1:t
    vel(i) = S(i+1)-S(i);
end

figure

plot(theta,vel)

end
