function I=keye(D,N)
% Identity operator.
%
% Input:
%   D = number of dimensions
%   N = scalar or vector (of size D), number of points in each dimension
%
% Output:
%   I = identity operator in separated format
    
    if length(N)==1
        N=N*ones(1,D);
    elseif length(N)~=D
        disp('idnt.m: N must have length 1 or length D.')
        I=[];
        return
    end
    
    U=cell(1,D);
    
    for d=1:D
        U{d}=reshape(speye(N(d)),[N(d)^2,1]);
    end
    
    I=ktensor(1,U);
    I=arrange(I);

end