function [ hMes ] = hRadar( x, rpos )
    % measurement function for the radar at rpos in 2D
    %  state x(6,1): x,vx,y,vy,theta,vtheta
    % rpos: (2,1)
    hMes = (x(1)-rpos(1))^2+(x(3)-rpos(2))^2;
end

