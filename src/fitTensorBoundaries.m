function [ tensorOutput, gridTnew, dx ] = fitTensorBoundaries( tens, gridT, meanT, covT, ngridT, lambda )
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


    % Check input
    validateattributes(tens,{'ktensor'},{'nonempty'})
    dim = ndims(tens);
    validateattributes(gridT,{'cell'},{'numel',dim})
    validateattributes(meanT,{'double'},{'numel',dim})
    validateattributes(ngridT,{'double'},{'numel',dim})
    validateattributes(covT,{'double'},{'size',[dim,dim]})
    validateattributes(lambda,{'double'},{'numel',1})
    
    
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

