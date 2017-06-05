function n = knumel(A)
% Number of elements in a ktensor

    n=0;
    for d=1:ndims(A)
        n = n + numel(A.U{d});
    end

end