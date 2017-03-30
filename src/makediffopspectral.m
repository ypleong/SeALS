function [n,grid,region,D,D2,fd1,fd2] = makediffopspectral(bdim,n,bcon,region)
% MAKEDIFFOP creates discretization grids and differentiation operators 
% in ktensor form and differentiation matrices.
% Inputs:
%   grid - the discretization grid.
%   n - number of grid points in each dimension.
%   bcon - the boundary conditions.
%   region - the boundaries of the goal region.
% Outputs:
%   n - number of grid points in each dimension.
%   grid - the discretization grid.
%   region - the boundaries of the goal region.
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

%% grid

size_bdim = size(bdim);
d = size_bdim(1);

if length(n(:)) == 1
    n = repmat(n,d,1);
end


%% region

if isempty(region) == 0
    
    % same boundary for all dimensions
    if length(region(:)) == 2
        %check if true dim
        if size(region) ~= [1 2]
            error('wrong size on region')
        end
        region = repmat(region,d,1);
    end   
end

%% create finite-difference matrices

fd1 = cell(1,d); %first order der
fd2 = cell(1,d); %second order der
grid = cell(1,d);

for i=1:d
    
    if bcon{i}{1} == 'p' %then bcon{i} = {'p'} so periodic
        [grid{i}, fd1{i}] = fourdif(n(i), 1);
        [grid{i}, fd2{i}] = fourdif(n(i), 2);
        
    else
        if ~isempty(region)
            [grid{i}, sd] = chebdifn(n(i), 2, bdim(i,1),bdim(i,2),region(i,1),region(i,2));
        else
            [grid{i}, sd] = chebdifn(n(i), 2, bdim(i,1),bdim(i,2));
        end
        fd1{i} = sd(:,:,1);
        fd2{i} = sd(:,:,2);
    end
    
end

if ~isempty(region)
    % point case for region, change region to closest point.
    if region(:,1) == region(:,2)
        for i=1:d
            gridi = grid{i};
            [~,index] = min(abs( gridi-region(i,1) ));
            region(i,:) = [gridi(index) gridi(index)];
        end
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

if isempty(region)
    
    for i=1:d
        DU = idUop;
        D2U = idUop;
        
        DU{i} = fd1{i}(:);
        D2U{i} = fd2{i}(:);
        
        D{i} = ktensor(DU);
        D2{i} = ktensor(D2U);
    end
    
else
    
    D = diffop_region_spectral(grid,fd1,region);
    D2 = diffop_region_spectral(grid,fd2,region);
    
end

end