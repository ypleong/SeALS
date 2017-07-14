function [bc] = makebcforward(bcon,n)
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
hasd = 0;

for i=1:d
    
    bconi = bcon{i};
    
    if strcmp(bconi{1},'p')
        bcU{i} = ones(n(i),1);
        
    elseif strcmp(bconi{1},'d')
        bcU{i} = create_boundary_dirichlet(bconi,n(i));
        hasd = 1;
        
    elseif strcmp(bconi{1},'n')
        bcU{i} = ones(n(i),1);
        
    elseif strcmp(bconi{1},'v') % vanishing boundary
        bcU{i} = ones(n(i),1);
        
    else
        error('wrong type of boundary condition')
    end
    
    
end

if hasd == 0
    bc = 0*oneTens(d,n);
else
    bc = ktensor(bcU);
end

bc = fixsigns(arrange(bc));

end

function [bc] = create_boundary_dirichlet(bcon,ndim)
    % create bc
    bc = zeros(ndim,1);
    bc([1 ndim],:) = bcon{2:3};
end