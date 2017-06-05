function one = kones(D,N)
% Generate a (dense) ones tensor.
% Input:
%   D = number of dimensions
%   N = scalar or D-vector specifying the size of each dimension

    U = cell(1,D);
    
    if length(N)==1
        N = N*ones(1,D);
    end

    for d=1:D
        U{d} = ones(N(d),1);
    end
    
    one = ktensor(1,U);
    one = normalize(one);

end