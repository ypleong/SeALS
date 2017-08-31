function [ tensorOutput, gridTnew ] = fitTensorBoundaries( tens, gridT, meanT, covT, ngridT, lambda )
%
%  Adjust the boundaries of the pdf to lambda sigmas around the mean
%
%       Input:
%  tens, 
%  gridT, 
%  meanT, 
%  covT, 
%  ngridT, 
%  lambda
%       Output:
%  tensorOutput, 
%  gridTnew 

    dim = ndims(tens);
    gridTnew = cell(dim,1);
    bdim = zeros(dim,2);
    dx = zeros(dim,1);
    for i=1:dim
        bdim(i,:) = [meanT(i)-lambda*sqrt(covT(i,i)) meanT(i)+lambda*sqrt(covT(i,i))];
        dx(i) = (bdim(i,2)-bdim(i,1))/(ngridT(i)-1);
        gridTnew{i} =  (bdim(i,1):dx(i):bdim(i,2))';
    end
    tensorOutput = interpTensor(tens,gridT,gridTnew);
    
end

