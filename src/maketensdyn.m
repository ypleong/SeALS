function [fTens,GTens,BTens,noise_covTens,qTens,RTens] = maketensdyn(f,G,B,noise_cov,q,R,x,grid)
% MAKETENSDYN creates ktensor representation of symbolic functions, given below.
% Inputs:
%   f - the function f in symbolic representation.
%   G - the function G in symbolic representation.
%   B - the function B in symbolic representation.
%   noise_cov - the noise covariance matrix.
%   q - the state cost in symbolic representation.
%   R - the control cost matrix.
%   x - the symbolic state vector.
%   grid - the discretization grid.
% Outputs: (matrices where the entries are ktensors)
%   fTens - f as ktensor.
%   GTens - G as ktensor.
%   BTens - B as ktensor.
%   noise_covTens - noise_cov as ktensor.
%   qTens - q as ktensor.
%   RTens - R as ktensor.
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 5 2015

%% create fTens
if isa(f,'sym') == 0 %not sym (constant matrix)
    f = sym(f); %convert to symbolic
end
fcell = fsym2fcell(f,x); %create cell representation
fTens = fcell2ftens(fcell,grid); %create ktensor

%% create GTens
if isa(G,'sym') == 0
    G = sym(G);
end
Gcell = fsym2fcell(G,x);
GTens = fcell2ftens(Gcell,grid);

%% create BTens
if isa(B,'sym') == 0
    B = sym(B);
end
Bcell = fsym2fcell(B,x);
BTens = fcell2ftens(Bcell,grid);

%% create noise_covTens
if isa(noise_cov,'sym') == 0
    noise_cov = sym(noise_cov);
end
noise_covcell = fsym2fcell(noise_cov,x);
noise_covTens = fcell2ftens(noise_covcell,grid);

%% create qTens
if isa(q,'sym') == 0
    q = sym(q);
end
qcell = fsym2fcell(q,x);
qTens = fcell2ftens(qcell,grid);

%% create RTens
if isa(R,'sym') == 0
    R = sym(R);
end
Rcell = fsym2fcell(R,x);
RTens = fcell2ftens(Rcell,grid);
 
end