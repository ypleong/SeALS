function [ ktens ] = ktensorGaussian( xMean, diagCov, gridT )
% Generate a ktensor from a diagonal gaussian pdf
%
% xMean  nx1
% diagCov nx1
% gridT cell(n,1)

    dim = length(xMean);
    validateattributes(gridT,{'cell'},{'numel',dim})
    validateattributes(diagCov,{'double'},{'numel',dim})
    
    p0vector = cell(dim,1);
    for i=1:dim
       p0vector{i} =  normpdf(gridT{i},xMean(i),sqrt(diagCov(i))); 
    end
    ktens= ktensor(p0vector);

end

