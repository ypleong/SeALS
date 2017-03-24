function [op] = makebcop(op,bcon,bsca,n,fd1)
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
%   bsca - boundary scaling. bsca(i,1) and bsca(i,2) is how much lower and
%   upper boundary in dimension i should be scaled.
%   n - the number of grid points in each dimension.
%   fd1 - cell of first order differentiation matrices. fd1{i} is
%   differentiation in dimension i.
% Output:
%   op - adjusted operator for boundary conditions as ktensor.

% Elis Stefansson, Aug 2015
% created for LHJB Toolbox

d = length(n);

for i=1:d
    
    bconi = bcon{i};
    
    if strcmp(bconi{1},'p')
        if bsca(i,1) ~= bsca(i,2)
            error('scaling for periodic boundary should be same for lower and upper boundary')
        end
        op = bc_periodic(op,i,bsca(i,:),n);
        
    elseif strcmp(bconi{1},'d')
        op = bc_dirichlet(op,i,bsca(i,:),n);
        
    elseif strcmp(bconi{1},'n')
        op = bc_neumann(op,i,bsca(i,:),n,fd1{i});
    else
        error('wrong type of boundary condition')
    end
    
end