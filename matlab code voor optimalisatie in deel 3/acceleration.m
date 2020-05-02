function acc = acceleration(S)

t = size(S);
theta = (1:1:t);
S = [S;S(1)];

vel = zeros(t);
for i = 1:t
    vel(i) = S(i+1)-S(i);
end
vel = [vel;vel(1)];

acc = zeros(t);
for i = 1:t
    acc(i) = vel(i+1)-vel(i);
end

% plot(theta,acc)

end
