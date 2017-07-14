function [bc] = makebcbackward(bc,bcon,n)
% MAKEBC modifies the ktensor to incooporate boundary conditions.
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
%   bc - adjusted ktensor for boundary conditions.

d = length(n);

for i=1:d
    
    bconi = bcon{i};
    
    if strcmp(bconi{1},'p')
        bc.U{i} = eliminate_boundary(bc.U{i},n(i));
        
    elseif strcmp(bconi{1},'d')
        continue
        
    elseif strcmp(bconi{1},'n')
        continue
        
    elseif strcmp(bconi{1},'v') % vanishing boundary
        continue
        
    else
        error('wrong type of boundary condition')
    end
    
    
end

bc = fixsigns(arrange(bc));

end

function [op] = eliminate_boundary(opU,ndim)
    % Eliminate whatever is happening on the boundary currently
    sr = size(opU,2);
    dimU = opU;
    dimU([1 ndim],:) = zeros(2,sr);
    op = dimU;
end