function plot2Dslice(F,slices_dim,coordinates,grid,lambda)
    
% Inputs:
% F - ktensor
% slices_dim - dimensions to plot
% coordinates - coordinates of other dimensions
% grid - cell of grid points

% Outputs:
% 3D surface plot

    if nargin < 5
        vv = 0;
    else
        vv = 1;
    end
    d = ndims(F);
    factors = F.lambda;
    
    for i = 1:d
        if all(i ~= slices_dim)
            factors = factors.*F.U{i}(coordinates(i),:)';
        end
    end
    
    kkksubU{1} = F.U{slices_dim(1)};
    kkksubU{2} = F.U{slices_dim(2)};
    if vv
        kkksub = -log(abs(double(ktensor(factors, kkksubU)))*lambda);
    else
        kkksub = abs(double(ktensor(factors, kkksubU)));
    end
    
    surf(grid{slices_dim(1)},grid{slices_dim(2)},kkksub','EdgeColor','none')
    xlabel(['x_{',num2str(slices_dim(1)),'}'])
    ylabel(['x_{',num2str(slices_dim(2)),'}'])
end