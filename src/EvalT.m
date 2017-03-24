function [y] = EvalT(tens,x,grid)
% EVALT evaluates a ktensor (function) at an individual point
% Inputs:
%   tens - the ktensor.
%   x - the point where the ktensor value should be evaluated.
%   grid - the discretization grid.
% Outputs:
%   y - the value at x.

d = length(x);

index = cell(d,1);
for i=1:d
    bas = grid{i};
    k = dsearchn(bas,x(i));
    index{i} = k;
end

y = tens(index{:});
end

