function [ y ] = EvalTMat( TMat, x, grid )
% EVALTMAT evaluates a matrix-valued ktensor function at a point x
% Inputs:
%   TMat - the matrix-valued ktensor function.
%   x - the points where the value of TMat should be evaluated.
%   grid - the discretization grid.
% Outputs:
%   y - the value at x (a matrix).

[n,m] = size(TMat);
y = zeros(n,m);

for i=1:n
    for j=1:m
        y(i,j) = EvalT(TMat{i,j},x,grid);
    end
end

end

