function [ meanT, covT ] = meanCovTensor( itens, gridT, weMean, weCov, weOnes )
% Compute the mean and covariance of a tensor
%
% itens, 
% gridT, cell(n,1)
% weMean, weights for the mean. cell(dim,dim);
% weCov, weights for the covariance, cell(dim,dim,dim);

    dim = ndims(itens);
    itens = itens*(1/intTens(itens, [], gridT, weOnes));
    meanT = zeros(dim,1);
    for i=1:dim
        meanT(i)  = intTens(itens, [], gridT, weMean(i,:));
    end
    covT = zeros(dim,dim);
    for i=1:dim
        for j = 1:dim
            covT(i,j) = intTens(itens, [], gridT, weCov(i,j,:)) - meanT(i)*meanT(j);
        end
    end

end

