function [P,P_mean] = power_cam

close all; clear;

%Excentricity = 0

    % Load our cam data
    cam1 = load('hefwet + geometrie zonder excentriciteit.mat');
    alpha = cam1.pressure_angle;
    R=sqrt(cam1.xpitch.^2+cam1.ypitch.^2)*0.001;
    omega = cam1.w;
    N = cam1.normalforce_tot;
    
    %Instantaneous power
    P_withoutexc = N.*sin(alpha).*R*omega;
    
    %Mean power
    P_withoutexc_mean = sum(P_withoutexc)/size(P_withoutexc,2);
    
%Excentricity =/= 0

    % Load our cam data
    cam2 = load('hefwet + geometrie.mat');
    alpha = cam2.pressure_angle;
    R=sqrt(cam2.xpitch.^2+cam2.ypitch.^2)*0.001;
    omega = cam2.w;
    e = cam2.exc*0.001;
    N = cam2.normalforce_tot;

    %Instantaneous power
    P_withexc=N.*(cos(alpha).*e+sin(alpha).*sqrt(R.^2-e^2))*omega;
    
    %Mean power    
    P_withexc_mean = sum(P_withexc)/size(P_withexc,2)
    
%Plot

    %With and without exc    
    figure
    subplot(1,2,1)
    plot(cam1.theta*180/pi,P_withoutexc);
    hold on
    plot(cam1.theta*180/pi,P_withoutexc_mean*ones(size(cam1.theta)));
    title("Power with centric follower")
    xlabel('\theta (°)');
    ylabel('Power (W)');
    subplot(1,2,2)
    plot(cam1.theta*180/pi,P_withexc);
    hold on
    plot(cam1.theta*180/pi,P_withexc_mean*ones(size(cam1.theta)));
    title("Power with excentric follower")
    xlabel('\theta (°)');
    ylabel('Power (W)');
    
    %Difference
    figure
    plot(cam1.theta*180/pi,P_withexc-P_withoutexc);
    title("Difference in powers")
    xlabel('\theta (°)');
    ylabel('Power (W)');
    
P = P_withexc;
P_mean = P_withexc_mean;

end