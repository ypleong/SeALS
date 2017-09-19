function [ outputCheck ] = checkGridFit( gridT, meanT, covT, lambdaMin, lambdaMax )
% Check if the given gaussian fits for the given grid

% gridT, meanT, covT, lambdaMin, lambdaMax

% It checks:
%   - the mean has to be withing lambdaMin times the cov to the border.
%   Otherwise the grid is too small. It checks for both borders
%   - the covariance cannot be smaller than the size of the grid divided by
%   lambdaMax. That avoids resolution problems with peaked pdfs.

    validateattributes(meanT,{'cell'},{'nonempty'})
    dim = length(gridT);
    validateattributes(meanT,{'double'},{'numel',dim})
    validateattributes(covT,{'double'},{'size',[dim,dim]})
    validateattributes(lambdaMin,{'double'},{'<', 10,'>', 0})
    validateattributes(lambdaMax,{'double'},{'<', 20,'>', 5})
        
    outputCheck = zeros(dim,1);
    for i=1:dim
       outputCheck(i) = any([meanT(i) - lambdaMin*sqrt(covT(i,i)) < gridT{i}(1) ...
                 ,  meanT(i) + lambdaMin*sqrt(covT(i,i)) > gridT{i}(end) ...
                 ,  (gridT{i}(end)- gridT{i}(1))/2 / sqrt(covT(i,i)) > lambdaMax]);         
    end

end

