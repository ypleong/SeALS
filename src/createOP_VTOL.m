function [ op ] = createOP_VTOL(D, D2, gridT, tol_err_op, u, tau, g, epsilon, x, dt, q )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    fFPE(1) =  x(2);
    fFPE(2) = -u*sin(x(5)) + epsilon*tau*cos(x(5));
    fFPE(3) =  x(4);
    fFPE(4) =  u*cos(x(5)) + epsilon*tau*sin(x(5)) - g;
    fFPE(5) =  x(6);
    fFPE(6) =  tau;
    
    explicit = 1;
    % Create Operator
    op = create_FP_op (  fFPE, q, dt, D, D2,gridT, tol_err_op, explicit,x);
end

