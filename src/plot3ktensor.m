function [  ] = plot3ktensor (xvector, inputTensor, dimOut, nOut )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


    ndim = size(inputTensor);
    figure
    hold on
    for i=1:nOut
       iOut = floor( i/nOut*ndim(3) );
       h_t(i) = pcolor (xvector{1},xvector{2}, inputTensor(:,:,iOut)');
       set(h_t(i), 'EdgeColor', 'none');
       set(h_t(i),'ZData',i*ones(ndim(2),ndim(1)))
       set(h_t(i),'FaceAlpha',0.1)
    end
    colorbar
    
    ndim = size(inputTensor);
    figure
    h_t = pcolor (xvector{1},xvector{2}, inputTensor(:,:,1)');
    set(h_t, 'EdgeColor', 'none');
    axis ([0,25,0,25,0,nOut])
    view([-37.5 -30])
    caxis([0 0.020])
    colorbar
    for i=1:nOut
       iOut = floor( i/nOut*ndim(3) );
       iOut
       set(h_t, 'CData' ,  inputTensor(:,:,iOut)' );
       
       set(h_t,'ZData',i*ones(ndim(2),ndim(1)))
       %set(h_t,'FaceAlpha',0.1)
       drawnow;
       pause(1.0/30.0);
    end
    

end

