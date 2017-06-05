function A = trncel(A,delta)
% Absolute truncation of factor matrices.  All entries less than delta in
% magnitude are removed.
%
% Input:
%   A = ctd
%   delta = truncation tolerance. If delta is empty or 0, do nothing.
%
% Output:
%   A = ctd

    % if delta empty or 0, do nothing
    if isempty(delta) || delta==0
        return
    end

    % truncate each dimension
    for d=1:ndims(A)
        A.U{d}=trncspmat(A.U{d},delta);
    end
    
    % re-normalize
    A = normalize(A);
    
end