function B = getterms(A,idx)
% extract the terms specificed by vector idx.  This is slow if you try to
% extract thousands of terms.  Ok for hundreds...
 
    B = cell(1,ndims(A));
    
    for d = 1:ndims(A)

        B{d} = A{d}(:,idx);
        
    end
    
    % extract weights
    B = ktensor(A.lambda(idx),B);

end