function [  ] = plot2DtimeMFibers( pk, timeVector, gridT )
% Plot the time evolution of the tensor marginalized for each dimension
%
% pk, 
% timeVector, %
% gridT,
%

    dim = ndims(pk{1});
    nT = arrayfun( @(x) length( gridT{x} ), 1:dim ); 
    [ weMean, weCov, weOnes ] = createWeights( gridT, nT );
    meanT = zeros(dim,length(timeVector));
    for k=1:length(timeVector)
        %pk{k} = pk{k}*(1/intTens(pk{k}, [], gridT, weOnes));
        [ meanT(:,k), ~ ] = meanCovTensor( pk{k}, gridT, weMean, weCov, weOnes );
    end
    figure
    title('Marginalized fibers for each dimension')
    for i=1:dim
        subplot(dim,1,i)
        fiberTime = zeros(length(gridT{i}),length(timeVector));
        for k=1:length(timeVector)
            tempT = intTens(pk{k}, [1:i-1 i+1:dim], gridT, []);
            fiberTime(:,k) = sum( repmat(tempT.lambda',size(tempT.U{i},1),1).*tempT.U{i},2);
        end
        handleOutput = pcolor(timeVector',gridT{i},fiberTime);
        shading interp 
        hold on
        plot(timeVector,meanT(i,:) )
        set(handleOutput,'EdgeColor','none');
        xlabel('time(s)')
        ylabel(['x_',num2str(i)])   
    end

end

