function [fFunc,GFunc,BFunc,noise_covFunc,qFunc,RFunc] = makefuncdyn(f,G,B,noise_cov,q,R,x)
%MAKEFUNCDYN creates MATLAB functions of symbolic functions, given below.
% Inputs:
%   f - the function f in symbolic representation.
%   G - the function G in symbolic representation.
%   B - the function B in symbolic representation.
%   noise_cov - the noise covariance matrix.
%   q - the state cost in symbolic representation.
%   R - the control cost matrix.
%   x - the symbolic state vector.
% Outputs:
%   fFunc - f as MATLAB function.
%   GFunc - G as MATLAB function.
%   BFunc - B as MATLAB function.
%   noise_covFunc - noise_cov as MATLAB function.
%   qFunc - q as MATLAB function.
%   RFunc - R as MATLAB function.
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 2015

%% create fFunc
if isa(f,'sym') == 0 %not sym (constant matrix)
    f = sym(f);
end
fFunc = matlabFunction(f,'Vars',{x});

%% create GFunc
if isa(G,'sym') == 0
    G = sym(G);
end
GFunc = matlabFunction(G,'Vars',{x});

%% create BFunc
if isa(B,'sym') == 0
    B = sym(B);
end
BFunc = matlabFunction(B,'Vars',{x});

%% create noise_covFunc
if isa(noise_cov,'sym') == 0
    noise_cov = sym(noise_cov);
end
noise_covFunc = matlabFunction(noise_cov,'Vars',{x});

%% create qFunc
if isa(q,'sym') == 0
    q = sym(q);
end
qFunc = matlabFunction(q,'Vars',{x});

%% create RFunc
if isa(R,'sym') == 0
    R = sym(R);
end
RFunc = matlabFunction(R,'Vars',{x});
 
end