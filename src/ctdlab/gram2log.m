function G = gram2log(A,B)
% gram matrix of tensors A and B.

    G = A.lambda * B.lambda';
    Gs = sign(G);
    G = log(abs(G) + eps);
    for d = 1:ndims(A)
      ABtmp = full(A{d}'*B{d});
      Gs = Gs .* sign(ABtmp);
      G = G + log(abs(ABtmp) + eps);
    end
    G = Gs .* exp(G);

end