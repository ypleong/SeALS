function [  ] = plot2d( xvector,zmatrix )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    hh = pcolor(xvector{1}',xvector{2}',zmatrix');
    colorbar
    set(hh, 'EdgeColor', 'none');
    xlabel('x')
    zlabel('p(x)')
    grid on


end

