function acc = acceleration(S)

t = size(S);
theta = (1:1:t);
S = [S;S(1)];

vel = velocity(S);
vel = [vel;vel(1)];

acc = zeros(t);
for i = 1:t
    acc(i) = vel(i+1)-vel(i);
end

% figure
% plot(theta,acc)

end
