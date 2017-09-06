function [ handleOutput ] = plot2DProjectionPoint(points, handleInput )

    sizeP = size(points);
    dim = sizeP(1);
    if nargin == 1
        handleOutput = {};        
        for i=1:dim 
            for j=(i+1):dim
                subplot(dim,dim,dim*(i-1)+j)
                hold on
                handleOutput{i,j} = scatter(points(i,:),points(j,:));
            end
        end    
    else
        for i=1:dim             
            for j=(i+1):dim
                set( handleInput{i,j},'XData', points(i,:),'Ydata',points(j,:));
            end
        end 
        drawnow limitrate
        pause(1.0/1000);        
    end

end