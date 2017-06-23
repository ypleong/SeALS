function idTens = oneTens(d, n)
% Create a d dimensional all ones tensor of size n

% Inputs:
% d - number of dimensions
% n - length of each dimension

% Output:
% idTens - ktensor of all ones

if length(n) == 1
    n = n*ones(1,d);
elseif length(n) ~= d
    error('length(n) must be the same as the dimension')
end

U = cell(d,1);
for i = 1:d
    U{i} = ones(n(i),1);
end

idTens = ktensor(U);

end