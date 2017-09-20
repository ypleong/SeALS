%% Create a ktensor with M terms in n dimensions
clear all
% common variables
dim = 3;
n = [101 121 131];
bdim = [0 10
        0 10
        0 10];
for i=1:dim
    dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1); % grid space
    gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))'; % grid vector
end
[ weMean, weCov, weOnes ] = createWeights( gridT, n );


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

%% boundary test

% Generate initial ktensor
meanT = [1 1 5 ]';
covT = diag([0.5 0.1 0.2]);
pcheck = ktensorGaussian( meanT, diag(covT), gridT );
plotkTensor(pcheck,gridT)

% Check boundary
lambdaMin = 4;
lambdaInitial = 6;
lambdaMax = 10;
checkGridFit(gridT, meanT, covT,lambdaMin,lambdaMax)

[ pcheckAfter, gridT_after, dx] = fitTensorBoundaries( pcheck, gridT, meanT, covT, n, lambdaInitial );
plotkTensor(pcheckAfter,gridT_after)


afigure
gridT = -6:0.2:6;
hold on
plot(gridT,normpdf(gridT,0,1))
plot(gridT,normpdf(gridT,0,6/4))
plot(gridT,normpdf(gridT,0,6/12),'.-')
legend('Original \lambda=6','Too Big \lambda=4','Too Small \lambda=12')


%% plotfibers test

% Show how the plotkTensorCell function works with 3 dimensions
clear all
dim = 3;
n = [101 121 131];
bdim = [0 10
        0 10
        0 10];
for i=1:dim
    dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1); % grid space
    gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))'; % grid vector
end

meanT = [1 1 5 ]';
covT = diag([0.5 0.1 0.2]);

t = 0:0.01:1;
pcheck = cell(length(t),1);
for k=1:length(t)
    pcheck{k} = (2+t(k))*ktensorGaussian( meanT+[t(k)*5,4*t(k)^2,sin(t(k)*2*pi)], (0.5+t(k))*diag(covT), gridT ) ... 
         +      ktensorGaussian( meanT+[4,7,5]'-[t(k)*5,3*t(k)^2,sin(t(k)*2*pi)], diag(covT), gridT );
    pcheck{k} = arrange(pcheck{k});
end
plotkTensorCell(pcheck,gridT,t,'plotSeq','marginalizedFiber')
plotkTensorCell(pcheck,gridT,t,'plotSeq','slider')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test als2 for diagonal gaussians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Common setup
clear all
dim = 2;
n = [101 101];
bdim = [-4 4
        -4 4];
    
for i=1:dim
    dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1); % grid space
    gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))'; % grid vector
end
meanT = [0 0]';

%% Default Gaussian
covT = diag([1 0.5^2]);
pcheck = ktensorGaussian( meanT , diag(covT), gridT );
plotkTensor(pcheck,gridT)
%% Test Angle
covT = diag([1 (1/2)^2]);
[X1,X2] = meshgrid(gridT{1}',gridT{2}');
theta = 0:0.02:pi/2;
nLambda = zeros(length(theta),1);
allLambda = cell(length(theta),1);
pk = cell(length(theta),1);
for i=1:length(theta)
    Rtheta = [cos(theta(i)) -sin(theta(i));sin(theta(i)) cos(theta(i))];
    gaussianXY = mvnpdf([X1(:) X2(:)],meanT',Rtheta*covT*Rtheta');
    gaussianXY = reshape(gaussianXY,length(gridT{2}),length(gridT{1}));
    [pKxy,Uo,out(i)] = cp_als(tensor(gaussianXY),101,'init','nvecs');
    [pKxy2, err, iter, e_list, t_step, illcond, noreduce] = als2(pKxy);
    pk{i} = pKxy2;
    %plotkTensor(pKxy2,gridT);
    allLambda{i} = pKxy2.lambda;
    nLambda(i) = length(pKxy2.lambda);
end
afigure
plot(theta*180/pi,nLambda)
xlabel('Rotation Angle')
ylabel('Number of terms')
grid on
hold on

figure
plotkTensorCell(pk,gridT,theta)

figure
plotkTensorCell(pk,gridT,theta,'plotSeq','marginalizedFiber')
% figure
% hold on
% for i=1:length(theta)
%    plot(allLambda{i},'Color',[i/length(theta) 0 1-i/length(theta)])
% end
% hold off


%% Just for the 45deg case, show one case
thetaI = pi/4;
Rtheta = [cos(thetaI) -sin(thetaI);sin(thetaI) cos(thetaI)];
gaussianXY = mvnpdf([X1(:) X2(:)],meanT',Rtheta*covT*Rtheta');
gaussianXY = reshape(gaussianXY,length(gridT{2}),length(gridT{1}));
[pKxy,Uo,out(i)] = cp_als(tensor(gaussianXY),101,'init','nvecs');
[pKxy2, err, iter, e_list, t_step, illcond, noreduce] = als2(pKxy,1e-4);
plotkTensor(pKxy2,gridT);

%% Check accuracy graph
epsV = logspace(-7,-1,200);
nLambda = zeros(length(epsV),1);
for i=1:length(epsV)
    [pKxy2, err, iter, e_list, t_step, illcond, noreduce] = als2(pKxy,epsV(i));
    nLambda(i) = length(pKxy2.lambda);
end
%afigure
semilogx(epsV,nLambda)
xlabel('Required Accuracy')
ylabel('Number of terms')
hold on
grid on
fitOut = fit(log10(epsV'),nLambda,'poly1')


%% Check covariance ratio
thetaI = pi/4;
Rtheta = [cos(thetaI) -sin(thetaI);sin(thetaI) cos(thetaI)];
rV = 0.1:0.001:1;
%rV = [2 4 8 16]; 
%rV = 1./rV;
nLambda = zeros(length(rV),1);
for i = 1:length(rV)
    covT = diag([1 rV(i)^2]);
    gaussianXY = mvnpdf([X1(:) X2(:)],meanT',Rtheta*covT*Rtheta');
    gaussianXY = reshape(gaussianXY,length(gridT{2}),length(gridT{1}));
    [pKxy,Uo,out(i)] = cp_als(tensor(gaussianXY),101,'init','nvecs');
    [pKxy2, err, iter, e_list, t_step, illcond, noreduce] = als2(pKxy);
    nLambda(i) = length(pKxy2.lambda);
end
afigure
plot(rV,nLambda)
xlabel('Covariance Ratio')
ylabel('Number of terms')
grid on


%% Check locus test
nLambda = zeros(100,1);
figure
hold on
for i=1:length(nLambda)
    [pKxy2, err, iter, e_list, t_step, illcond, noreduce] = als2(pKxy,1e-8);
    nLambda(i) = length(pKxy2.lambda);
    plotkTensor(pKxy2,gridT);
end
xlabel('Number of terms')
grid on

nLambda = length(pKxy2.lambda);
filterVector = zeros(nLambda,1);
for i=1:nLambda
    figure(i)
    
    filterVector(i) = 1;
    plotkTensor(filterBasis(pKxy2,filterVector),gridT);
    
end


%% 3D gaussian

% Common setup
clear all
dim = 3;
n = [101 101 101];
bdim = [-4 4 
        -4 4
        -4 4];
    
for i=1:dim
    dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1); % grid space
    gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))'; % grid vector
end
meanT = [0 0 0]';

%% Default Gaussian
covT = diag([1 0.5^2 0.25^2]);
pcheck = ktensorGaussian( meanT , diag(covT), gridT );
plotkTensor(pcheck,gridT)
%% Test Angle
[X1,X2,X3] = meshgrid(gridT{1}',gridT{2}',gridT{3}');
omega = 0:0.1:pi/2;
phi = 0:0.1:pi/2;
nLambda = zeros(length(omega),length(phi));
for i=1:length(omega)
    for j=1:length(phi)
        Rtheta = euler2rot( omega(i), phi(j), 0 );
        %Rtheta = [cos(theta(i)) -sin(theta(i));sin(theta(i)) cos(theta(i))];
        gaussianXY = mvnpdf([X1(:) X2(:) X2(:)],meanT',Rtheta*covT*Rtheta');
        gaussianXY = reshape(gaussianXY,length(gridT{2}),length(gridT{1}),length(gridT{3}));
        [pKxy,Uo,out(i)] = cp_als(tensor(gaussianXY),101,'init','nvecs');
        [pKxy2, err, iter, e_list, t_step, illcond, noreduce] = als2(pKxy);
        %pk{i} = pKxy2;
        %plotkTensor(pKxy2,gridT);
        nLambda(i,j) = length(pKxy2.lambda);
    end
end
figure
surf(omega*180/pi,phi*180/pi,nLambda)
xlabel('\Omega')
ylabel('\phi')
zlabel('Number of terms')
grid on
hold on

