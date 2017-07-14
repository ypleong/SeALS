function [op] = makebcopforward(op,bcon,n)
% MAKEBCOP modifies the operator to incooporate boundary conditions.
% Input:
%   bcon - the boundary conditions. bcon{i} is boundary conditions in
%   dimension i:
%      bcon{i} = {'p'}. Periodic
%      bcon{i} = {'d',val_lo,val_up}. Dirichlet with val_lo and val_up 
%      for lower and upper boundary in dimension i.
%      bcon{i} = {'n',val_lo,val_up}. Neumann with val_lo and val_up for
%      lower and upper boundary in dimension i.
%      (where all the values are for the desirability function.)
%   n - the number of grid points in each dimension.
% Output:
%   op - adjusted operator for boundary conditions as ktensor.

d = length(n);

for i=1:d
    
    bconi = bcon{i};
    
    if strcmp(bconi{1},'p')
        continue
        
    elseif strcmp(bconi{1},'d')
        op.U{i} = eliminate_boundary(op.U{i},n(i));
        
    elseif strcmp(bconi{1},'n')
        continue
        
    elseif strcmp(bconi{1},'v') % vanishing boundary
        op.U{i} = eliminate_boundary(op.U{i},n(i));
        
    else
        error('wrong type of boundary condition')
    end
    
    
end

op = fixsigns(arrange(op));

end

function [op] = eliminate_boundary(opU,ndim)
    % Eliminate whatever is happening on the boundary currently
    sr = size(opU,2);
    dimU = opU;
    for jj = 1:sr
        el = reshape(dimU(:,jj),ndim,ndim);
        el([1 ndim],:) = zeros(2,ndim);
        dimU(:,jj) = el(:);
    end
    op = dimU;
end