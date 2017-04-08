function [op,bc] = incorpregion(op,bc,region,grid,regval,regsca)
% INCORPREGION modifies operator and boundary conditions to incooporate
% Dirichlet boundary conditions at the goal region.
% Inputs:
%   op - the operator as ktensor
%   bc - the boundary conditions as ktensor.
%   region - the dimensions of the goal region.
%   grid - the discretization grid.
%   regval - the Dirichlet value at the goal region.
%   regsca - the boundary scaling at the goal region.
% Outputs:
%   op - the adjusted operator.
%   bc - the adjusted boundary conditions.

% Elis Stefansson, Aug 2015

sr = ncomponents(op);
sr_b = ncomponents(bc);
d = length(grid);

% select out the region from the operator
opAdd = cell(1,d);
opRem = op;

bcAdd = cell(d,1); %will be the added bc in region.
bcd = bc; %will be a copy of bc in region

for i=1:d
    
    gridi = grid{i};
    ni = length(gridi);
    d_ind = find( (gridi >= region(i,1)) & (gridi <= region(i,2)));
    d_out = find( (gridi < region(i,1)) | (gridi > region(i,2)));
    
    %Go through each rank and eliminate nodes outside region
    temp = reshape(opRem.U{i},ni,ni*sr);
    temp(d_out,:) = 0;
    opRem.U{i} = reshape(temp,ni*ni,sr);
    
    %Create identity for nodes inside region
    id = eye(ni);
    id(d_out,:) = 0;
    id(:,d_out) = 0;
    opAdd{i} = id(:);
    
    %Go through the boundary and eliminate nodes outside region
    newBCU = bcd.U{i};
    newBCU(d_out,:) = zeros(length(d_out),sr_b);
    bcd.U{i} = newBCU;
    bcAdd{i} = zeros(ni,1);
    bcAdd{i}(d_ind,1) = ones(length(d_ind),1);
    
end

newbcop = regsca*ktensor(opAdd);

op = op - opRem + newbcop;

bc = bc - bcd + regsca*regval*ktensor(bcAdd);

end