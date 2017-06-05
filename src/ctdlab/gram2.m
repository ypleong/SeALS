function G = gram2(A,B)
% gram matrix of tensors A and B.

    G = A.lambda * B.lambda';
    for d = 1:ndims(A)
        G = G .* (A{d}'*B{d});
    end

end