function [y] = EvalT(tens,x,grid)
% EVALT evaluates a ktensor (function) at an individual point
% Inputs:
%   tens - the ktensor.
%   x - the point where the ktensor value should be evaluated.
%   grid - the discretization grid.
% Outputs:
%   y - the value at x.

y = tens.lambda';
for i=1:length(x);
    y = y.*interp1( grid{i}, tens.U{i}, x(i) );
end

y = sum(y);
end

