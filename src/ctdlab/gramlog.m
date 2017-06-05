function G = gramlog(A)
% Gram matrix of a tensor A.  Normally you might exploit symmetry, but
% because of sparse data structure it's not clear that this will be
% faster...
    
    G = A.lambda * A.lambda';   
    Gs = sign(G);
    G = log(abs(G) + eps);
    for d = 1:ndims(A)
      Atmp = A{d}'*A{d};
      Gs = Gs .* sign(Atmp);
      G = G + log(abs(Atmp) + eps);
    end
    G = Gs .* exp(G);
end