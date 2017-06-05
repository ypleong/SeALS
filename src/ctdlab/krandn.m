function B=krandn(D,N,R,density)
% Create random sparse ctd with normally distributed elements.
%
% Input:
%   D = number of dimensions
%   N = scalar or D-vector of dimension sizes
%   R = separtion rank
%   density = scalar between 0 and 1, approximate fraction of non-zeros.
%
% Output:
%   B = a random sparse ktensor

    % if N is a scalar, make all dimensions of this size
    if length(N)==1
        N=N*ones(1,D);
    end

    B=cell(1,D);
    for d=1:D
        if density==1
            B{d}=randn(N(d),R);
        else
            B{d}=sprandn(N(d),R,density);
        end
    end
    B=ktensor(ones(R,1),B);
    
    % normalize and sort
    B = arrange(B);

end