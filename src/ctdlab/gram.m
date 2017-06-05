function G = gram(A)
% Gram matrix of a tensor A.  Normally you might exploit symmetry, but
% because of sparse data structure it's not clear that this will be
% faster...
    
    G = A.lambda * A.lambda';    
    for d = 1:ndims(A)
        G = G .* (A{d}'*A{d});
    end

end