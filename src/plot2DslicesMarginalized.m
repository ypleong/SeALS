function [ handleOutput ] = plot2DslicesMarginalized(Ftensor, point, gridT, handleInput )
    % Plots a visualization of the Tensor. For each pair of dimensions it
    % plots a 2D slice around the selected point. The value of the
    % dimensions that are not sliced is kept constant the the vlue of the
    % selected point
    %
    % Inputs:
    % F - ktensor (n dims)
    % point - (nx1 vector) 
    % grid - cell of grid points of {n vectors}
    % handleInput - handles to generate animation
    
    % Outputs:
    % Slice Representation
    %
    % Usage: 
    %  - Static Plot: call the function without the handleInput
    %  - Animation: load the plot first using a static plot and then call
    %  the function in a loop, example code:
    %      
    %      figure
    %      handleSlices = plot2DslicesAroundPoint( Ftensor{1}, point(1,:), gridT);
    %      for k = 2:30:length(t)
    %          plot2DslicesAroundPoint( Ftensor{k}, point(:,k), gridT, handleSlices)    
    %      end
    %
    %
    
    dim = ndims(Ftensor);
    indexT = zeros(dim,1);
    for i=1:dim
       [~, indexT(i)] = min( abs(point(i)-gridT{i}));
    end
    if nargin == 3
        handleOutput = {};
        
        for i=1:dim 
            handleOutput{i,i} = subplot(dim,dim,dim*(i-1)+i);
            grid on
            hold on
            xlabel(['x_',num2str(i)])
            title('     Basis Function')
            sizeU = size(Ftensor.U{i});
            for k = 1:sizeU(2)
                plot(gridT{i},Ftensor.lambda(k)*Ftensor.U{i}(:,k),'black');
            end
            for j=(i+1):dim
                subplot(dim,dim,dim*(i-1)+j)
                handleOutput{i,j} = plot2DsliceM(Ftensor,[i,j],gridT);
            end
        end    
    else
        for i=1:dim             
            handle_r = subplot(dim,dim,dim*(i-1)+i);
            cla(handle_r)
            %ylim([-1,1])
            hold on
            sizeU = size(Ftensor.U{i});
            for k = 1:sizeU(2)
                plot(gridT{i},Ftensor.lambda(k)*Ftensor.U{i}(:,k),'b');
                
            end
            for j=(i+1):dim
                plot2DsliceM(Ftensor,[i,j],indexT,gridT,handleInput{i,j});
            end
        end 
        drawnow limitrate
        pause(1.0/1000);
        
    end
        
        


end