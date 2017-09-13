function [ handleoutput ] = plotkTensor( F, gridT, varargin )
% Plot a ktensor using different plotting styles
%
% F: ktensor with 'dim' dimensions
% gridT: a dim,1 cell with the grid coordinates for each dimension
% varargin: input arguments to change the default behaviour

validateattributes(F,{'ktensor'},{'nonempty'})
dim = ndims(F);
validateattributes(gridT,{'cell'},{'size',[1,dim]})
ngrid = arrayfun(@(x) length(gridT{x}),1:dim);
[weMean,weCov,weOnes] = createWeights(gridT,ngrid);
[ meanT, ~ ] = meanCovTensor( F, gridT, weMean, weCov, weOnes );

params = inputParser;

%addRequired(params,'shape',@(x)validateattributes(x,charchk,nempty))
addParameter(params,'takeLog',0,@isscalar);


defaultPlot2D = 'slide';
expectedPlot2D = {'no','marginalized','slide'};
addParameter(params,'plot2D',defaultPlot2D,...
                 @(x) any(validatestring(x,expectedPlot2D)));
             
defaultPlot1D = 'marginalized';
expectedPlot1D = {'no','marginalized','basis'};
addParameter(params,'plot1D',defaultPlot1D,...
                 @(x) any(validatestring(x,expectedPlot1D)));
             
             
defaultplotType = 'surf';
expectedplotType = {'surf','pcolor'};
addParameter(params,'plotType',defaultplotType,...
                 @(x) any(validatestring(x,expectedplotType)));
             
addParameter(params,'slidePoint',meanT,@(x) length(x)==dim)
addParameter(params,'handleInput',[],@(x) iscell(x) && size(x)==[dim,dim])



params.parse(varargin{:});

%% Copy from params object
plot2D = params.Results.plot2D;
plot1D = params.Results.plot1D;
plotType = params.Results.plotType;
takeLog = params.Results.takeLog;
slidePoint = params.Results.slidePoint;
handleInput = params.Results.handleInput;

%% Plot ktensor
switch plot2D
    case 'no'
        
    case 'slide'
        [ handleOutput ] = plot2DslicesAroundPoint(F, slidePoint, gridT, handleInput, plotType);
    case 'marginalized'
        [ handleOutput ] = plot2DslicesMarginalized(F, gridT, handleInput);
end

end

