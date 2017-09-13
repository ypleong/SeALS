function [handleOutput] = plot2Dslice(F,slices_dim,coordinates,gridT, handleInput,lambda,plotType)
    
% Inputs:
% F - ktensor (n dims)
% slices_dim - (2x1 vector) dimensions to plot
% coordinates - (nx1 vector) coordinates of other dimensions
% grid - cell of grid points

% Outputs:
% 3D surface plot

    d = ndims(F);
    factors = F.lambda;
    
    for i = 1:d
        if all(i ~= slices_dim)
            factors = factors.*F.U{i}(coordinates(i),:)';
        end
    end
    
    kkksubU{1} = F.U{slices_dim(1)};
    kkksubU{2} = F.U{slices_dim(2)};
    if ~isempty(lambda)
        kkksub = -log(abs(double(ktensor(factors, kkksubU)))*lambda);
    else
        kkksub = (double(ktensor(factors, kkksubU)));
    end
    
    if isempty(handleInput)
        switch plotType
            case 'surf'
                handleOutput = surf(gridT{slices_dim(1)},gridT{slices_dim(2)},kkksub','EdgeColor','none');
            case 'pcolor'
                handleOutput = pcolor(gridT{slices_dim(1)},gridT{slices_dim(2)},kkksub');
                set(handleOutput,'EdgeColor','none');
        end
        view(0,90)
        %caxis([0 0.3])
		xlabel(['x_{',num2str(slices_dim(1)),'}'])
		ylabel(['x_{',num2str(slices_dim(2)),'}'])
    else
        switch plotType
            case 'surf'
                set( handleInput, 'ZData', kkksub' )
            case 'pcolor'
                set( handleInput, 'CData', kkksub' )
        end
    end
        
end