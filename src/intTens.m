function [Fint] = intTens(F, w, dim)

% F - ktensor
% w - integration weight in each dimension in cell
% dim - dimensions to integrate over

nd = ndims(F);
if length(dim) ~= nd
    Fintc = cell(1, nd - length(dim));
else
    Fintc = cell(1, 1);
end

Fint_lambda = F.lambda';
int_all = 0;

for ii = 1:nd
    if any(ii == dim)
        Fintc{ii} = w{ii}*F.U{ii};
        Fint_lambda = Fint_lambda.*Fintc{ii};
    else
        Fintc{ii} = F.U{ii}; 
    end
    int_all = int_all + size(Fintc{ii},1);
end

if int_all == nd
    Fint = sum(Fint_lambda);
else
    Fint = ktensor(F.lambda, Fintc);
end