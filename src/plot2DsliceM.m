function [handleOutput] = plot2DsliceM(itens,slices_dim,gridT, handleInput,lambda)
%
%   Plot coordinates slices_dim and marginalize all others
%
% Inputs:
% F - ktensor (n dims)
% slices_dim - (2x1 vector) dimensions to plot
% grid - cell of grid points

% Outputs:
% 3D surface plot
    dim = ndims(itens);
    i = slices_dim(1);
    j = slices_dim(2);
    if (dim>2)
        tempT = intTens(itens, [1:i-1 i+1:j-1 j+1:dim], gridT, []);
    else
        tempT = itens;
    end
    ktempT = ktensor(tempT.lambda,tempT.U{i},tempT.U{j});
    if nargin == 3
        handleOutput = pcolor(gridT{i},gridT{j},double(ktempT)');
        set(handleOutput,'EdgeColor','none');
		xlabel(['x_{',num2str(slices_dim(1)),'}'])
		ylabel(['x_{',num2str(slices_dim(2)),'}'])
    else 
       set( handleInput, 'CData', double(ktempT)' )
    end
        
end