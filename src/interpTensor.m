function [ tensOutput ] = interpTensor( tens, oldgridT, newgridT )
%
% Change the grid of a tensor and interpolate it to the new grid
% Extrapolate using zeros.
%
% - tens - the ktensor.
% -oldgridT, the old grid as a cell with size dim of vectors
% -newgridT, the new grid as a cell with size dim of vectors, they might 
% have different size than the oldgrid
    
    validateattributes(tens,{'ktensor'},{'nonempty'})
    dim = ndims(tens);
    validateattributes(gridT,{'cell'},{'numel',dim})
    validateattributes(newgridT,{'cell'},{'numel',dim})
    
    tensOutput = tens; % copy tensor
    kterms = length(tens.lambda);
    
    for i=1:dim
        for j=1:kterms
            tensOutput.U{i}(:,j) = interp1(oldgridT{i},tens.U{i}(:,j),newgridT{i},'PCHIP',0);
        end
    end

end

