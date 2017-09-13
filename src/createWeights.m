function [ weMean, weCov, weOnes ] = createWeights( gridT, n )
    % Create weights for integration
    dim = length(gridT);
    weMean = cell(dim,dim);
    for i=1:dim
        for j = 1:dim
            if i == j
                weMean{i,j} = gridT{j};
            else
                weMean{i,j} = ones(n(j),1);
            end         
        end
    end
    weOnes = cell(dim,1);
    for i=1:dim
        weOnes{i} = ones(n(i),1);
    end
    weCov = cell(dim,dim,dim);
    for i=1:dim
        for j = 1:dim
            weCov(i,j,:) = weOnes;
            if i == j
                weCov{i,j,i} = gridT{j}.*gridT{j};
            else
                weCov{i,j,i} = gridT{i};
                weCov{i,j,j} = gridT{j};
            end         
        end
    end

end

