function [D,D2,fd1,fd2] = makediffop(grid,n,h,ord_of_acc,bcon,region)
% MAKEDIFFOP creates differentiation operators in ktensor form and
% differentiation matrices.
% Inputs:
%   grid - the discretization grid.
%   n - number of grid points in each dimension.
%   h - the grid step in each dimension.
%   ord_of_acc - oder of accuracy for the derivatives.
%   bcon - the boundary conditions.
%   region - the dimensions of the goal region.
% Outputs:
%   D - cell of first order differentiation operators. D{i} is
%   differentiation in dimension i.
%   D2 - cell of second order differentiation operators. D2{i} is
%   differentiation in dimension i.
%   fd1 - cell of first order differentiation matrices. fd1{i} is
%   differentiation in dimension i.
%   fd2 - cell of second order differentiation matrices. fd2{i} is
%   differnetiation in dimension i.
%
%   See also MAIN_RUN

% Elis Stefansson, Aug 5 2015

%% create finite-difference matrices

d = length(grid);

fd1 = cell(1,d); %first order der
fd2 = cell(1,d); %second order der

if isempty(region) == 0 %have region
    fd1_edg = cell(1,d); %edge derivatives for region
    fd2_edg = cell(1,d);
end

for i=1:d
    
    if length(bcon{i}) == 1 %then bcon{i} = {'p'} so periodic
        
        fd1{i} = fdmatrix(n(i),1,ord_of_acc(1,:),1);
        fd1{i} = 1/h(i)*fd1{i};
        fd2{i} = fdmatrix(n(i),2,ord_of_acc(2,:),1);
        fd2{i} = 1/(h(i)^2)*fd2{i};
        
    else
        
        fd1{i} = fdmatrix(n(i),1,ord_of_acc(1,:));
        fd1{i} = 1/h(i)*fd1{i};
        fd2{i} = fdmatrix(n(i),2,ord_of_acc(2,:));
        fd2{i} = 1/(h(i)^2)*fd2{i};
        
    end
    
    if isempty(region) == 0
        
        fd1_edg{i} = fdmatrix(n(i),1,ord_of_acc(1,:));
        fd1_edg{i} = 1/h(i)*fd1_edg{i};
        fd2_edg{i} = fdmatrix(n(i),2,ord_of_acc(2,:));
        fd2_edg{i} = 1/(h(i)^2)*fd2_edg{i};
        
    end
    
end

%% create differentiation operators

% create identity operator
idUop = cell(1,d);
for i=1:d
    id = eye(n(i));
    idUop{i} = id(:);
end

D = cell(1,d);
D2 = cell(1,d);

if isempty(region) == 1
    
    for i=1:d
        DU = idUop;
        D2U = idUop;
        
        fd1i = fd1{i};
        fd2i = fd2{i};
        
        DU{i} = fd1i(:);
        D2U{i} = fd2i(:);
        
        D{i} = ktensor(DU);
        D2{i} = ktensor(D2U);
    end
    
else
    
    D = diffop_region(grid,fd1,region,ord_of_acc,fd1_edg,1);
    D2 = diffop_region(grid,fd2,region,ord_of_acc,fd2_edg,2);
    
end

end