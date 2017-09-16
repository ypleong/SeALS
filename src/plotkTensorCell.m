function [ handleoutput ] = plotkTensorCell( Fcell, gridT, seqVector, varargin )
% Plot a sequence of ktensors

params = inputParser;
% Check required parameters
validateattributes(Fcell,{'cell'},{'nonempty'},'ktensor');
validateattributes(Fcell{1},{'ktensor'},{'nonempty'})

dim = ndims(Fcell{1});
validateattributes(gridT,{'cell'},{'numel',dim})
if (isempty(seqVector))
    seqVector = 1:length(Fcell);
else
    validateattributes(seqVector,{'double'},{'numel',length(Fcell)})
end
    
% Check optional parameters
defaultPlotSeq = 'slider';
expectedPlotSeq = {'marginalizedFiber','video','slider','snapshots'};
addParameter(params,'plotSeq',defaultPlotSeq,...
                 @(x) any(validatestring(x,expectedPlotSeq)));
addParameter(params,'skip',10,@(x) isnumeric(x)&& x>0 && x<length(Fcell))
addParameter(params,'sliderArrow',0.02,@(x) isscalar(x) && x>0 && x<1)

params.StructExpand = false;
opt = struct;
opt.fileName = 'saveDefault.avi';
opt.saveVideo = 0;

addParameter(params,'options',opt,@(x) isstruct(x))
addParameter(params,'nSnapshots',8,@(x) isscalar(x) && isinteger(x) && x>1 && x<length(Fcell))
%Parse parameters
params.parse(varargin{:});
plotSeq = params.Results.plotSeq;
skipFrames = params.Results.skip;
options = params.Results.options;
arrowJumpSlider = params.Results.sliderArrow;
nSnapshots = params.Results.nSnapshots;
%% Plot sequence

switch plotSeq
    case 'marginalizedFiber'
        plot2DtimeMFibers( Fcell, seqVector, gridT )
    case 'snapshots'
        for k=1:nSnapshots
            figure(k)
            indexT = max(round(k*length(seqVector)/nSnapshots),1);
            handleSlices = plot2DslicesMarginalized( Fcell{indexT}, gridT,[]);
        end
    case 'video'
        hf = figure;
        if (options.saveVideo)
                title('Probability Density Function Evolution');
                writerObj = VideoWriter('pdf_gaussian_.avi');
                writerObj.Quality = 100;
                %set(hf,'Visible','on');
                open(writerObj);
                set(gcf,'Renderer','OpenGL'); %to save to file
        end
        handleSlices = plot2DslicesMarginalized( Fcell{1}, gridT,[]);
        for k = 2:skipFrames:length(seqVector)
            plot2DslicesMarginalized( Fcell{k}, gridT, handleSlices);
            pause(1.0/100);
            if (options.saveVideo)   
                writeVideo(writerObj, getframe(gcf));
            end  
        end
        if (options.saveVideo)   
            close(writerObj);
        end  
    case 'slider'
        h_figure = figure;
        handleSlices = plot2DslicesMarginalized( Fcell{1}, gridT,[]);
        b_text = uicontrol('Parent',h_figure,'Style','text','Position',[51,50,30,23]);
        callbackSliderAA = @(source,event) callbackSlider(source,event,Fcell,gridT,handleSlices,b_text);
        
        b_slider = uicontrol('Parent',h_figure,'Style','slider','Position',[81,54,419,23],...
                'value',1, 'min',1, 'max',length(seqVector), 'SliderStep',[max(1/(length(seqVector)-1),arrowJumpSlider) max(1/(length(seqVector)-1),0.1)],'Callback', callbackSliderAA);
end


    function callbackSlider(source,~,Fcell,gridT,handleSlices,b_text)
        sliderValue = round(source.Value);
        source.Value = sliderValue;
        set(b_text,'String',num2str(sliderValue));
        plot2DslicesMarginalized( Fcell{sliderValue}, gridT, handleSlices);
    end
    
end