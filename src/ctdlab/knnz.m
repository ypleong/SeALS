function n=knnz(A)
% Count the non-zero elements of a spktensor
        
    n=0;
    for d=1:ndims(A)
        n=n+nnz(A.U{d});
    end

end