function [ tensOutput ] = interpTensor( tens, oldgridT, newgridT )
% Change the grid of a tensor
%  tens - the ktensor.
    tensOutput = tens; % copy tensor
    kterms = length(tens.lambda);
    dim = ndims(tens);
    for i=1:dim
        for j=1:kterms
            tensOutput.U{i}(:,j) = interp1(oldgridT{i},tens.U{i}(:,j),newgridT{i},'PCHIP',0);
        end
    end

end

