function [ xdot ] = fVTOL( x, u, tau, epsilon, g )
% x:[6,1] state ( x, xdot, y, ydot, theta, thetadot)
% u and tau are control inputs
% epsilon is a parameter, ~0.01
% g is gravity 9.8
    xdot = zeros(6,1);

    xdot(1) =  x(2);
    xdot(2) = -u*sin(x(5)) + epsilon*tau*cos(x(5));
    xdot(3) =  x(4);
    xdot(4) =  u*cos(x(5)) + epsilon*tau*sin(x(5)) - g;
    xdot(5) =  x(6);
    xdot(6) =  tau;


end

