function [ handleOutput ] = plot2DslicesAroundPointGaussian(covPK, point, gridT, handleInput )
    % Plots a visualization of the Tensor. For each pair of dimensions it
    % plots a 2D slice around the selected point. The value of the
    % dimensions that are not sliced is kept constant the the vlue of the
    % selected point
    %
    % Inputs:
    % covPK - covariance matrix (nxn vector)
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
    %      handleSlices = plot2DslicesAroundPoint( Phat(:,:,1), point(1,:), gridT);
    %      for k = 2:30:length(t)
    %          plot2DslicesAroundPoint( Phat(:,:,k), point(:,k), gridT, handleSlices)    
    %      end
    %
    %
    
    dim = length(point);
    if nargin == 3
        handleOutput = {};
        for i=1:dim 
            for j=(i+1):dim
                subplot(dim,dim,dim*(i-1)+j)
                matMin = zeros(dim,2);
                matMin(i,1) = 1;
                matMin(j,2) = 1;
                covPKmin = matMin'*covPK*matMin;
                [X1,X2] = meshgrid(gridT{i},gridT{j});
                handleOutput{i,j} = pcolor(gridT{i},gridT{j},reshape( ...
                    mvnpdf([X1(:) X2(:)], [point(i),point(j)],covPKmin),length(gridT{j}),length(gridT{i})));
                h = colorbar;
                set(handleOutput{i,j}, 'EdgeColor', 'none');
                xlabel(['x_',num2str(i)])
                ylabel(['x_',num2str(j)])
                axis tight
                grid on
            end
        end    
    else
        for i=1:dim             
            for j=(i+1):dim
                matMin = zeros(dim,2);
                matMin(i,1) = 1;
                matMin(j,2) = 1;
                covPKmin = matMin'*covPK*matMin;
                [X1,X2] = meshgrid(gridT{i},gridT{j});
                set(handleInput{i,j},'XData',gridT{i}, ...
                                     'YData',gridT{j}, ...
                                     'CData',reshape(  ...
                    mvnpdf([X1(:) X2(:)], [point(i),point(j)],covPKmin),length(gridT{j}),length(gridT{i})));
                                     % 'XLim', [gridT{i}(1), gridT{i}(end)], ...
                                     % 'YLim', [gridT{j}(1), gridT{j}(end)]);
                %axis(handleInput{i,j},'auto')
                %axis([gridT{i}(1), gridT{i}(end), gridT{j}(1), gridT{j}(end)]);
            end
        end 
        drawnow limitrate
        pause(1.0/1000);
        
    end
        
        


end

