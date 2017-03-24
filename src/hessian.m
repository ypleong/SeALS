function [hess] = hessian(D, D2)
% HESSIAN generates the Hessian operator as a cell array with ktensor
% elements where the ktensors corresponds to the entires in the usual
% Hessian matrix.
% Input: 
%   D - cell with first order derivatives operators. D{i} is
%   differentiation with respect to dimension i.
%   D2 - cell with first order derivatives operators. D2{i} is
%   differentiation with respect to dimension i.
% Output:
%   hess - the Hessian operator.
%
% See also MAKEOP

d = ndims(D{1});
hess = cell(d,d);

for i=1:d
    for j=i:d
        
        if i==j
            hess{i,j} = D2{i};
        else
            hess{i,j} = SRMultM(D{i},D{j});
        end
        
        if i ~= j
            hess{j,i} = hess{i,j};
        end
        
    end
end

end

