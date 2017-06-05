function b = mvm(A,x)
% Matrix-vector multiplication

    D = ndims(A);
    N2 = size(A);
    N = sqrt(N2);
    rA = length(A.lambda);
    rx = length(x.lambda);
    
    b.U = cell(1,D);
    b.lambda = zeros(rA*rx,1);
    
    for d = 1:D
        % i need to change this to use sparse indexing....
        b.U{d} = spalloc(N(d),rA*rx,N(d)*(rA*rx));
        ii = 1;
        for r = 1:rA
            Atmp = reshape(A.U{d}(:,r),[N(d),N(d)]);
            for s = 1:rx
                b.U{d}(:,ii) = Atmp * x.U{d}(:,s);
                b.lambda(ii) = A.lambda(r)*x.lambda(s);
                ii = ii+1;
            end
        end
    end
    
    b = ktensor(b.lambda,b.U);
    b = arrange(b);

end