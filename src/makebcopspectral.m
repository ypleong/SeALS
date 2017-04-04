function [op] = makebcopspectral(op,bcon,bsca,n,fd1)
% MAKEBCOP modifies the operator to incooporate boundary conditions.
% Input:
%   bcon - the boundary conditions. bcon{i} is boundary conditions in
%   dimension i:
%      bcon{i} = {'p'}. Periodic. Do nothing.
%      bcon{i} = {'d',val_lo,val_up}. Dirichlet with val_lo and val_up 
%      for lower and upper boundary in dimension i.
%      bcon{i} = {'n',val_lo,val_up}. Neumann with val_lo and val_up for
%      lower and upper boundary in dimension i.
%      (where all the values are for the desirability function.)
%   bsca - boundary scaling. bsca(i,1) and bsca(i,2) is how much lower and
%   upper boundary in dimension i should be scaled.
%   n - the number of grid points in each dimension.
%   fd1 - cell of first order differentiation matrices. fd1{i} is
%   differentiation in dimension i.
% Output:
%   op - adjusted operator for boundary conditions as ktensor.


d = length(n);

for i=1:d
    
    bconi = bcon{i};
    
    if strcmp(bconi{1},'p')
        continue
    elseif strcmp(bconi{1},'d')
        op = bc_dirichlet(op,i,bsca(i,:),n);
        
    elseif strcmp(bconi{1},'n')
        op = bc_neumann(op,i,bsca(i,:),n,fd1{i});
        
    elseif strcmp(bconi{1},'v') % vanishing boundary
        % Eliminate whatever is happening on the boundary currently
        sr = ncomponents(op);
        dimU = op.U{i};
        ndim = n(i);
        for jj=1:sr
            el = reshape(dimU(:,jj),ndim,ndim);
            el(1,:) = zeros(1,ndim);
            %el(ndim,:) = zeros(1,ndim); %not needed
            dimU(:,jj) = el(:);
        end
        op.U{i} = dimU;
        op = arrange(op);
    else
        error('wrong type of boundary condition')
    end
    
end