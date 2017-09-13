%% Create a ktensor with M terms in n dimensions
clear all
% common variables
dim = 3;
n = [101 101 101];
bdim = [0 10
        0 10
        0 10];
for i=1:dim
    dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1); % grid space
    gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))'; % grid vector
end
[ weMean, weCov, weOnes ] = createWeights( gridT, n )


x0 = [2 2 2]';
diagSigma = [0.3 0.2 0.2];
p0 = ktensorGaussian( x0, diagSigma, gridT );
[ meanT, covT ] = meanCovTensor( p0, gridT, weMean, weCov, weOnes )
plot2DslicesAroundPoint( p0, x0, gridT, [],'pcolor');

x1 = [8 7 8]';
diagSigma1 = [0.3 0.1 0.2];
p1 = ktensorGaussian( x1, diagSigma1, gridT );
[ meanT, covT ] = meanCovTensor( p1, gridT, weMean, weCov, weOnes )
plot2DslicesAroundPoint( p1, x1, gridT, [],'pcolor');

pplus = p0+p1+ktensorGaussian( [5 3 5 ], [0.5 0.1 0.2], gridT );
[ meanT, covT ] = meanCovTensor( pplus, gridT, weMean, weCov, weOnes )
plot2DslicesAroundPoint( pplus, meanT, gridT, [],'pcolor');

pplus1 = arrange(pplus*(1/intTens(pplus, [], gridT, weOnes)));
[ meanT, covT ] = meanCovTensor( pplus1, gridT, weMean, weCov, weOnes )



figure
title('Marginalized fibers for each dimension')
Uu = cell(dim,1);
for i=1:dim
    subplot(dim,1,i)
    Uu{i} = sum( repmat(pplus1.lambda',size(pplus1.U{i},1),1).*pplus1.U{i},2);  
    plot(gridT{i},Uu{i}) 
    grid on
    xlabel(['x_',num2str(i)])
    ylabel(['U_',num2str(i)])   
end
pplus3 = ktensor(Uu);

% 
afigure
for i=1:dim
    subplot(dim,1,i)
    tempT = intTens(pplus1, [1:i-1 i+1:dim], gridT, weOnes);
    tempTT = zeros(length(gridT{1}),1);
    for k=1:length(pplus1.lambda)
        tempTT = tempTT +pplus1.lambda(k)*tempT{i}(:,k);
    end
    plot(gridT{i},tempTT)
end


figure
for i=1:dim    
    for j=(i+1):dim
        subplot(dim,dim,dim*(i-1)+j);
        tempT = intTens(pplus1, [1:i-1 i+1:j-1 j+1:dim], gridT, weOnes);
        handleOutput = pcolor(gridT{i},gridT{j},double(ktensor(tempT.lambda,tempT.U{i},tempT.U{j}))');
        set(handleOutput,'EdgeColor','none');
    end
end

