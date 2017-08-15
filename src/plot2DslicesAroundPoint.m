function [ handleOutput ] = plot2DslicesAroundPoint(Ftensor, point, gridT, handleInput )
    % Plots a visualization of the Tensor. For each pair of dimensions it
    % plots a 2D slice around the selected point. The value of the
    % dimensions that are not sliced is kept constant the the vlue of the
    % selected point
    %
    % Inputs:
    % F - ktensor (n dims)
    % point - (nx1 vector) 
    % grid - cell of grid points of {n vectors}

    % Outputs:
    % Slice Representation

    dim = ndims(Ftensor);
    indexT = zeros(dim,1);
    for i=1:dim
       [~, indexT(i)] = min( abs(point(i)-gridT{i}));
    end
    if nargin == 3
        handleOutput = {};
        %figure
        for i=1:dim 
            handleOutput{i,i} = subplot(dim,dim,dim*(i-1)+i);
            hold on
            xlabel(['x_',num2str(i)])
            title('Basis Function')
            sizeU = size(Ftensor.U{i});
            for k = 1:sizeU(2)
                plot(gridT{i},Ftensor.U{i}(:,k));
            end
            for j=(i+1):dim
                subplot(dim,dim,dim*(i-1)+j)
                handleOutput{i,j} = plot2Dslice(Ftensor,[i,j],indexT,gridT);
            end
        end    
    else
        for i=1:dim             
            handle_r = subplot(dim,dim,dim*(i-1)+i);
            cla(handle_r)
            ylim([-1,1])
            hold on
            sizeU = size(Ftensor.U{i});
            for k = 1:sizeU(2)
                plot(gridT{i},Ftensor.U{i}(:,k));
                %set(handleInput{i,i,k}, 'YData', Ftensor.U{i}(:,k) );
            end
            for j=(i+1):dim
                %subplot(dim,dim,dim*(i-1)+j)
                plot2Dslice(Ftensor,[i,j],indexT,gridT,handleInput{i,j});
            end
        end 
        drawnow limitrate
        pause(1.0/1000);
        
    end
        
        


end

