function [  ] = plot3ktensor (xvector, inputTensor, printOption, nOut  )
% Display a 3Way tensor
%
% xvector: {nx1, mx1} 2cell of vectors for x and y
% inputTensor size=(n,m,p)  3 variables
% printOption: string to determine how to print the slices
%              'static': print all the slices
%              'points': scatter3
%              'slider': change the slice with a ui slider
%              'time': change the slices in a video 
% nOut: number of slices. recomended: 10 for static and 100 for dynamic
    
    if ndims(inputTensor) ~= 3 
        disp(['InputTensor is not a 3way tensor but a ',num2str(ndims(inputTensor)),'way tensor'])
        return
    end
    
    if nargin < 3 || isempty(printOption)
       printOption = 'slider' 
    end

    if nargin < 4 || isempty(nOut)
        switch printOption
            case 'static'
                nOut = 10;
            case 'slider'
                nOut = 100;
            case 'time'
                nOut = 100;
        end
    end

    ndim = size(inputTensor);
    f = figure;
    hold on
    switch printOption
        case 'static'

            for i=1:nOut
               iOut = floor( i/nOut*ndim(3) );
               h_t(i) = pcolor (xvector{1},xvector{2}, inputTensor(:,:,iOut)');
               set(h_t(i), 'EdgeColor', 'none');
               set(h_t(i), 'ZData',i*ones(ndim(2),ndim(1)))
               set(h_t(i), 'FaceAlpha',0.1)
            end
            view([-37.5 30])
            axis ([xvector{1}(1),xvector{1}(end),xvector{2}(1),xvector{2}(end),xvector{1}(1),xvector{3}(end)])
            caxis([0 max(max(max(inputTensor)))])
            colorbar

        case 'slider'
            h_t = pcolor (xvector{1},xvector{2}, inputTensor(:,:,1)');
            set(h_t, 'EdgeColor', 'none');
            axis ([xvector{1}(1),xvector{1}(end),xvector{2}(1),xvector{2}(end),xvector{1}(1),xvector{3}(end)])
            view([-37.5 +30])
            caxis([0 max(max(max(inputTensor)))])
            colorbar

            callbackSliderAA = @(source,event) callbackSlider(source,event,h_t,ndim);

            b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
                  'value',ndim(3)/2, 'min',0, 'max',nOut, 'Callback', callbackSliderAA);
        case 'time'    
            h_t = pcolor (xvector{1},xvector{2}, inputTensor(:,:,1)');
            set(h_t, 'EdgeColor', 'none');
            axis ([xvector{1}(1),xvector{1}(end),xvector{2}(1),xvector{2}(end),xvector{1}(1),xvector{3}(end)])
           
            view([-37.5 +30])
            caxis([0 max(max(max(inputTensor)))])
            colorbar
            for i=1:nOut
               iOut = floor( i/nOut*ndim(3) );
               set(h_t, 'CData' ,  inputTensor(:,:,iOut)' );
               set(h_t, 'ZData',i*ones(ndim(2),ndim(1)))
               drawnow;
               pause(1.0/30.0);
            end
        case 'points'
            maxValue = max(max(max(inputTensor)));
            counter = 0;
            for i=1:length(xvector{1})
                for j=1:length(xvector{2})
                    for k=1:length(xvector{3})
                        if (inputTensor(i,j,k)>maxValue/10)
                           counter = counter + 1;
                           xyzC(counter,:) = [xvector{1}(i),xvector{2}(j),xvector{3}(k),inputTensor(i,j,k)];                           
                        end
                    end
                end
            end
            h_t = scatter3(xyzC(:,1),xyzC(:,2),xyzC(:,3),0.2,xyzC(:,4));
            axis ([xvector{1}(1),xvector{1}(end),xvector{2}(1),xvector{2}(end),xvector{1}(1),xvector{3}(end)])
           
            view([-37.5 +30])
            grid on
            axis equal
    end
     
    function callbackSlider(source,event,h_t,ndim)
        val = source.Value;
        iOut = floor( val/nOut*ndim(3) );
        set(h_t, 'CData' ,  inputTensor(:,:,iOut)' );
        set(h_t, 'ZData',xvector{3}(iOut)*ones(ndim(2),ndim(1)))
    end
        
        
    

end

