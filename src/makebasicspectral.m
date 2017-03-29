function [noise_cov,R] = makebasicspectral(noise_cov,R,G,B)
% MAKEBASIC initalize basic variables. For example, if ord_of_acc is
% specified as a scalar (same order of accuracy in all dimensions) then
% MAKEBASIC converts ord_of_acc into the general form. 
% Inputs:
%   noise_cov - the variance of the noise.
%   R - the control cost matrix (may be scalar).
%   G - the function G in symbolic representation.
%   B - the function B in symbolic representation.
%   See MAIN_PROGRAM for a detailed description.
% Outputs:
%   noise_cov - the noise covariance as a matrix.
%   R - the control cost matrix as a matrix.
%
% See also, MAIN_RUN.


%% noise covariance
sizeB = size(B);
if length(noise_cov) == 1
    noise_cov = noise_cov*eye(sizeB(2));
else %then vector for the diagonal elements. 
    noise_cov = diag(noise_cov);
end

%% matrix R
sizeG = size(G);
if length(R) == 1
    R = R*eye(sizeG(2));
end
end