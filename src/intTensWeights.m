function [w] = intTensWeights(n, bdim, gridt)

% n - number of grid points in each dimension. 
% bdim - the dimensions of the hyperrectangle domain; bdim(i,1) and
%   bdim(i,2) is the lower and upper boundary in dimension i.
% gridt - grid type 'p' - periodic (fourier), 'c' - chebyshev

% NOTE: DOES NOT INCLUDE INTERIOR BC

nd = length(n);

w = cell(1,nd);

for ii = 1:nd
    if gridt(ii) == 'p'
        w{ii} = 2*pi*ones(1,n(ii))/n(ii);
    elseif gridt(ii) == 'c'
        w{ii} = clencurt(n(ii), bdim(ii,1), bdim(ii,2));
    end
end
        