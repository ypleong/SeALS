function [ tensOutput ] = interpTensor( tens, oldgridT, newgridT )
%
% Change the grid of a tensor and interpolate it to the new grid
% Extrapolate using zeros.
%
% - tens - the ktensor.
% -oldgridT, the old grid as a cell with size dim of vectors
% -newgridT, the new grid as a cell with size dim of vectors, they might 
% have different size than the oldgrid

    tensOutput = tens; % copy tensor
    kterms = length(tens.lambda);
    dim = ndims(tens);
    for i=1:dim
        for j=1:kterms
            tensOutput.U{i}(:,j) = interp1(oldgridT{i},tens.U{i}(:,j),newgridT{i},'PCHIP',0);
        end
    end

end

