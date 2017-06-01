function plot2Dslice(F,slices_dim,coordinates,grid)
    
% Inputs:
% F - ktensor
% slices_dim - dimensions to plot
% coordinates - coordinates of other dimensions
% grid - cell of grid points

% Outputs:
% 3D surface plot

    
    d = ndims(F);
    factors = Fsub.lambda;
    
    for i = 1:d
        if all(i ~= slices_dim)
            factors = factors.*Fsub.U{i}(coordinates(i),:)';
        end
    end
    
    kkksubU{1} = Fsub.U{slices_dim(1)};
    kkksubU{2} = Fsub.U{slices_dim(2)};
    kkksub = double(ktensor(factors, kkksubU));
    
    figure; 
    surf(grid{slices_dim(1)},grid{slices_dim(2)},kkksub,'EdgeColor','none')
    xlabel(['x_',num2str(slices_dim(1))])
    ylabel(['x_',num2str(slices_dim(2))])
end