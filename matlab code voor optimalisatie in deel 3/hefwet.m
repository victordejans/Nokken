function S = hefwet
%returns a vector with the values of our hefwet (as defined in hefwet.mat) per degree


%init
theta = (1:1:360)';
S = zeros(size(theta));


%defining hefwet for each degree
for i = 1:60
    S(i) = 20*(theta(i)/60-1/pi*sin(pi*theta(i)/60)); %cycloid1
end
for i = 61:105
    S(i) = 20 + 15*(theta(i-60)/45+1/pi*sin(pi*theta(i-60)/45)); %cycloid2
end
for i = 106:150 
    S(i) = 35; %horizontal
end
for i = 151:290
    S(i) = 35*(1-theta(i-150)/140+1/(2*pi)*sin(2*pi*theta(i-150)/140)); %cycloid6
end


%some values calculated by this matlab code are smaller than zero with a very small absolute value. these are of course incorrect. we correct them to zero
S = abs(S);


% plot(theta,S)


end