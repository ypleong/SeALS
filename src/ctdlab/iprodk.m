function x = iprodk(A,B,k)
% Inner product of tensor A of arbitrary rank, and ktensor B of rank 1, 
% excluding dimension k.  This function is used for fixed point iteration.
% Input:
%   A = ktensor whose rank-1 approximation is sought
%   B = ktensor of rank 1
%   k = dimension to collapse around
% Output:
%   x = vector of size size(A,k)

    % check that B is of rank 1
    if length(B.lambda) > 1
        fprintf('iprodk:  B must be rank 1.\n')
        x=[];
        return
    end

    % iterate through all dimensions except k
    dims = 1:ndims(A);
    dims(k) = [];
    
    % collapse
    %coefMat = A.lambda;
    coefMat = A.lambda';
    for d = dims
        %coefMat = coefMat .* (A{d}'*B{d});
        tmp = full(B{d}' * A{d});
        coefMat = coefMat .* tmp;
    end
    coefMat = coefMat';
    
    % take linear combination of vectors in dimension k
    x = A{k} * coefMat;

end