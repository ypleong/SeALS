function A = poswts(A)
% make weights all positive.

    rnk = length(A.lambda);

    ineg = find(A.lambda<0);
    
    I = ones(rnk,1);
    I(ineg) = -1;
    
    I = spdiags(I,0,rnk,rnk);
    
    A{1} = A{1} * I;
    
    A.lambda(ineg) = -A.lambda(ineg);


end