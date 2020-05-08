function power

close all; clear;

%Excentricity = 0

    % Load our cam data
    cam = load('hefwet + geometrie zonder excentriciteit.mat');

    alpha = cam.pressure_angle;
    R=sqrt(cam.xpitch.^2+cam.ypitch.^2)*0.001;
    omega = cam.w;
    N = cam.normalforce_tot;

    P_withoutexc = N.*sin(alpha).*R*omega;
    
    figure
    plot(cam.theta*180/pi,P_withoutexc);
    title("Power with centric follower")
    xlabel('\theta (°)');
    ylabel('Power (W)');
    
%Excentricity =/= 0

    % Load our cam data
    cam2 = load('hefwet + geometrie.mat');

    alpha = cam2.pressure_angle;
    R=sqrt(cam2.xpitch.^2+cam2.ypitch.^2)*0.001;
    omega = cam2.w;
    e = cam2.exc*0.001;
    N = cam2.normalforce_tot;

    P_withexc=N.*(cos(alpha).*e+sin(alpha).*sqrt(R.^2-e^2))*omega;
    
    
    figure
    subplot(1,2,1)
    plot(cam.theta*180/pi,P_withoutexc);
    title("Power with centric follower")
    xlabel('\theta (°)');
    ylabel('Power (W)');
    subplot(1,2,2)
    plot(cam.theta*180/pi,P_withexc);
    title("Power with excentric follower")
    xlabel('\theta (°)');
    ylabel('Power (W)');
    
%Difference

    figure
    plot(cam.theta*180/pi,P_withexc-P_withoutexc);
    title("Difference in powers")
    xlabel('\theta (°)');
    ylabel('Power (W)');

end