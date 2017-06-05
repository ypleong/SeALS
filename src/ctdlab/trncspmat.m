function A = trncspmat(A,delta)
% For a sparse matrix, truncate all elements less than delta using a mask.  
% There might be a better way to do this, but this way avoids having 
% explicit indexing.
% Input:
%   A = sparse matrix
%   delta = cutoff threshold
% Ouput:
%   A = truncated input

    [i,j,s] = find(abs(A)>delta);
    [m,n] = size(A);
    mask = sparse(i,j,s,m,n);
    A = mask .* A;
    
end