function [ handleoutput ] = plotkTensorCell( Fcell, gridT, seqVector, varargin )
% Plot a sequence of ktensors


validateattributes(Fcell,{'cell'},{'nonempty'});
validateattributes(Fcell{1},{'ktensor'},{'nonempty'})

dim = ndims(Fcell{1});
validateattributes(gridT,{'cell'},{'size',[1,dim]})
if (isempty(seqVector))
    seqVector = 1:length(Fcell);
end

params = inputParser;

defaultPlotSeq = 'marginalizedFiber';
expectedPlotSeq = {'marginalizedFiber','video','slider'};
addParameter(params,'plotSeq',defaultPlotSeq,...
                 @(x) any(validatestring(x,expectedPlotSeq)));
addParameter(params,'skip',10,@(x) isnumeric(x)&& x>0 && x<length(Fcell))
params.parse(varargin{:});
plotSeq = params.Results.plotSeq;
skipFrames = params.Results.skip;

%% Plot sequence

switch plotSeq
    case 'marginalizedFiber'
        plot2DtimeMFibers( Fcell, seqVector, gridT )
    case 'video'
        handleSlices = plot2DslicesMarginalized( Fcell{1}, gridT,[]);
        for k = 2:skipFrames:length(seqVector)
            plot2DslicesMarginalized( Fcell{k}, gridT, handleSlices);
            pause(1.0/100);
        end
    case 'slider'
        h_figure = figure;
        handleSlices = plot2DslicesMarginalized( Fcell{1}, gridT,[]);
        b_text = uicontrol('Parent',h_figure,'Style','text','Position',[51,50,30,23]);
        callbackSliderAA = @(source,event) callbackSlider(source,event,Fcell,gridT,handleSlices,b_text);
        
        b_slider = uicontrol('Parent',h_figure,'Style','slider','Position',[81,54,419,23],...
                'value',1, 'min',1, 'max',length(seqVector), 'SliderStep',[max(1/(length(seqVector)-1),0.02) max(1/(length(seqVector)-1),0.1)],'Callback', callbackSliderAA);
end


    function callbackSlider(source,~,Fcell,gridT,handleSlices,b_text)
        sliderValue = round(source.Value);
        source.Value = sliderValue;
        set(b_text,'String',num2str(sliderValue));
        plot2DslicesMarginalized( Fcell{sliderValue}, gridT, handleSlices);
    end
    
end