function Ffunc = external_forces

%init
theta = (1:1:360)';
Ffunc = zeros(size(theta));


for i = 1:60
    Ffunc(i) = 0;
end

for i = 61:110
    Ffunc(i) = 150/50*(i-60); 
end

for i = 111:160
    Ffunc(i) = 250;
end

for i = 161:250
    Ffunc(i) = -230;
end

for i = 251:360
    Ffunc(i) = 0;
end

% plot(theta,Ffunc)

end