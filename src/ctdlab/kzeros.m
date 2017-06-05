function z = kzeros(D,N)
% Generate a zero sparse ktensor.
% Input:
%   D = number of dimensions
%   N = scalar or D-vector of dimension sizes

    if length(N)==1
        N = N*ones(1,D);
    end
    
    U = cell(1,D);
    for d=1:D
        U{d} = zeros(N(d),1);
    end
    
    z = ktensor(0,U);

end