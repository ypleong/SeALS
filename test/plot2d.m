function [  ] = plot2d( xvector,zmatrix )
% Plot a 2d mesh using pcolor.
%
% xvector: {nx1,mx1} cell with the x and y vectors
% zmatrix nxm matrix with values

    if ndims(zmatrix) == 2
        hh = pcolor(xvector{1}',xvector{2}',zmatrix');
        colorbar
        set(hh, 'EdgeColor', 'none');
        xlabel('x')
        ylabel('')
        grid on        
    else
        disp('Zmatrix is not a 2way matrix')
    end




end

