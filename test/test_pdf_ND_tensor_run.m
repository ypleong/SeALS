
%% Test Tensor Form PDF Evolution ND
clear all

%% Initialization, User-defined variables

caseStr = '2d_pos';
    
if strcmp('ND_forward',caseStr) 
     % This is just a constant velocity dynamics in N dimensions
     
    dim = 5; % change it!
    x = sym('x',[dim,1]); %do not change
    fprintf ('Running case  "%s"  with %d dimensions\n', caseStr, dim)
   
    measure = 0;
    fFPE = zeros(dim,1);
    for i=1:dim
        fFPE(i) = i/2;
        n(i) = 101+10*i;
        x0(i) = i;
        bdim(i,:) = [-2 8*i];
        bcon{i} = {'d',0,0};
        diagSigma(i) = 0.1*i;
        qdiag(i) = 0.1;
    end
    xhat = x0;


elseif strcmp('pendulum2D',caseStr) 
    
    % Linear pendulum
    
    dim = 2;
    x = sym('x',[dim,1]); %do not change
    fprintf ('Running case  "%s"  with %d dimensions\n', caseStr, dim)
    measure = 0;
    fFPE = [x(2);-x(1)];
    n = [201 201];
    x0 = [0.2, 0.01];
    diagSigma = [0.03 0.001];
    bdim = [-pi pi
           -10 10];
    bcon = { {'p'}, {'d',0,0}};
    
elseif strcmp('pendulum2D_NL',caseStr) 
     % This is the classical non-linear pendulum in 2d
     
    dim = 2;
    x = sym('x',[dim,1]); %do not change
    fprintf ('Running case  "%s"  with %d dimensions\n', caseStr, dim)
   
    measure = 0;
    fFPE = [x(2);-sin(x(1))];
    n = [301 301];
    x0 = [0.2, 0.01];
    diagSigma = [0.03 0.001];
    bdim = [-pi pi
           -10 10];
    bcon = { {'p'}, {'d',0,0}};
    
elseif strcmp('2d_pos',caseStr) 
     % This case is the simple double integrator (constant speed)
    % It checks the use of smaller measurements than the size of the state
    % vector
    
    dim = 2;
    x = sym('x',[dim,1]); %do not change
    fprintf ('Running case  "%s"  with %d dimensions\n', caseStr, dim)       
    measure = 1;
    mesureFreq = 1000;
    HK = [1 0];
    
    fFPE = [x(2);0];
    n = [201 301];
    x0 = [0.0, 0.50];
    x0hat = [0.0, 0.00];
    diagSigma = [0.05 0.4];
    bdim = [-4 5
           -4 4];
    bcon = { {'d',0,0}, {'d',0,0}};
    qdiag = [0.001 0.1];
    
elseif strcmp('3d_pos',caseStr)
    % This case is the simple integrator in the third dimension     
    
    dim = 3;
    x = sym('x',[dim,1]); %do not change
    fprintf ('Running case  "%s"  with %d dimensions\n', caseStr, dim)
    
    measure = 1;
    mesureFreq = 1000;
    HK = [1 0 0];
    fFPE = [x(2);x(3);0];
    n = [201 201 301];
    x0 = [0.0, 0.0, 0.50]';
    x0hat = [0.0, 0.0, 0.00]';
    diagSigma = [0.1 0.2 0.3];
    bdim = [-4 6
           -4 4
           -4 4];
    bcon = { {'d',0,0}, {'d',0,0},{'d',0,0}};
    qdiag = [0.001 0.05 0.02];
    
elseif strcmp('nd_hypersphere',caseStr) 
    dim = 3;
    x = sym('x',[dim,1]); %do not change
    fprintf ('Running case  "%s"  with %d dimensions\n', caseStr, dim)
    % This case shows the ability of the tensor framework to recover the
    % projection of a hypershere on a plane of the same dimensions
    
    measure = 0;
    mesureFreq = 1000;
    HK = [1 0 0];
    
    fFPE = [x(2);x(3);0];
    n = [201 201 301];
    x0 = [0.0, 0.0, 0.50]';
    x0hat = [0.0, 0.0, 0.00]';
    diagSigma = [0.1 0.2 0.3];
    bdim = [-4 6
           -4 4
           -4 4];
    bcon = { {'d',0,0}, {'d',0,0},{'d',0,0}};
    qdiag = [0.001 0.05 0.02];
elseif (strcmp('nonLinear4D',caseStr))
    dim = 4;
    x = sym('x',[dim,1]); %do not change   
    fFPE = [x(2);x(3);x(4);-x(3)*(kpen1+kpen2)-sin(x(1))*kpen1*kpen2];
    kpen1 = (2*pi/3)^2; 
    measure = 0;
else
    error('No case was specified')
end

% Consistency checks
if (length(fFPE) ~= dim )
   error('the dimensions of the dynamics are not "dim"') %check
end
% Derived Variables
if (~exist('x0hat','var'))
   x0hat = x0; 
end
for i=1:dim
    dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1); % grid space
    gridT{i} =  (bdim(i,1):dx(i):bdim(i,2))'; % grid vector
end
sigma02Matrix = diag(diagSigma); % [0.8 0.2; 0.2 sigma02];
for i=1:dim
   p0vector{i} =  normpdf(gridT{i},x0hat(i),sqrt(diagSigma(i))); 
end


%% Model parameters

if (~exist('qdiag','var') )
    qdiag = zeros(dim,1);
    for i=1:dim
        qdiag(i) = 0.3; 
    end
end
q =  diag( qdiag );

if (~exist('mesureFreq','var'))
   mesureFreq = 1; 
end


% Simulation Parameters
dt = 0.001;
finalt = 2;
t = 0:dt:finalt;

pk = cell(length(t),1);
zkcompressed = cell(length(t),1);
pk{1} =  ktensor(p0vector);
xGT = zeros(dim,length(t));
xkalman = zeros(dim,length(t));
covKalman = zeros(dim,dim,length(t));
expec = zeros(dim,length(t));
cov = zeros(dim,dim,length(t));
pz = cell(dim,1);

% Step up variables
cov(:,:,1) = sigma02Matrix;
expec(:,1) = x0hat;
% visualize, check boundaries
plot2DslicesAroundPoint( pk{1}, expec(:,1), gridT);

xGT(:,1) = x0;
xkalman(:,1) = x0hat;
covKalman(:,:,1) = sigma02Matrix;


%% Tensor parameters
if (~exist('bcon','var') )
    for i=1:dim
        bcon{i} = {'d',0,0};
    end
end
bsca = []; %no manual scaling
region = [];
regval = 1;
regsca = [];
sca_ver = 1;

tol_err_op = 1e-5;
tol_err = 1e-6;
als_options = [];
als_variant = []; %{10,50};
debugging = 0;
explicit = 1;

fFPE_function = matlabFunction(fFPE,'Vars',x);
fFPEdiff = jacobian(fFPE,x);
fFPEdiff_function = matlabFunction(fFPEdiff,'Vars',x);

for i=1:dim
    acc(:,i) = [2,2]';
end
[D,D2,~,~] = makediffop(gridT,n,dx,acc,bcon,region);
op = create_FP_op ( fFPE, q, dt, D, D2,gridT, tol_err_op, explicit,x);
[ weMean, weCov, weOnes ] = createWeights( gridT, n );

%% Iterate over 

time1 = tic;
for k = 2:length(t)
    
    if ( explicit == 1 )
        pk{k} = SRMultV( op, pk{k-1});
        if length(pk{k}.lambda) > 1
            [pk{k},~] = tenid(pk{k},tol_err_op,1,9,'frob',[],fnorm(pk{k}),0);

            pk{k} = fixsigns(arrange(pk{k}));
            [pk{k}, err_op, iter_op, enrich_op, t_step_op, cond_op, noreduce] = als2(pk{k},tol_err_op);

        end
    else 
        [pk{k}, err, iter, Fcond, enrich, t_step, illcondmat, maxit, maxrank, F_cell, B_cell, b_cell] = ...
        als_sys(op,pk{k-1},[],tol_err,als_options,debugging);
    end
    
    % Get Mean and Covariance values
    for i=1:dim
        expec(i,k)  = intTens(pk{k}, [], gridT, weMean(i,:));
    end
    for i=1:dim
        for j = 1:dim
            cov(i,j,k) = intTens(pk{k}, [], gridT, weCov(i,j,:)) - expec(i,k)*expec(j,k);
        end
    end

    % Integrate the GT
    xGT(:,k)  = deval(ode45( @(t,x) tempCall(fFPE_function,t,x),[t(k-1) t(k)], xGT(:,k-1)'),t(k));
    
    % Integrate KF
    xkalman(:,k)  = deval(ode45( @(t,x) tempCall(fFPE_function,t,x),[t(k-1) t(k)], xkalman(:,k-1)'),t(k));
    stm = (eye(dim) + fFPEdiff_function(xkalman(:,k-1))*dt);
    covKalman(:,:,k) = stm*covKalman(:,:,k-1)*stm'+dt*q;
    

    %% Measure
    if (measure && mod(k,mesureFreq)==0 )
        Rmeas = diag(HK*diag(diagSigma)*HK');
        mMes = size(HK);
        zmes = HK*xGT(:,k) + randn(mMes(1),1).*sqrt(Rmeas)';

        for i=1:dim
            for j=1:length(zmes)
                if ( HK(j,i) == 1 )
                    pz{i} =  normpdf(gridT{i},zmes(j),sqrt(Rmeas(j,j))); 
                elseif (all(HK(:,i)==0) )
                    pz{i} = ones(n(i),1);
                end
            end
        end

        zkcompressed{k} = ktensor(pz);

        % Direct PDF Bayesian Measurement
        pk{k} = HadTensProd(pk{k},zkcompressed{k});
        pk{k} = pk{k} *(1/  intTens(pk{k}, [], gridT, weOnes));
        
        for i=1:dim
            expecAf(i,k)  = intTens(pk{k}, [], gridT, weMean(i,:));
        end
        for i=1:dim
            for j = 1:dim
                covAf(i,j,k) = intTens(pk{k}, [], gridT, weCov(i,j,:)) - expecAf(i,k)*expecAf(j,k);
            end
        end

        % Kalman Measurement
        SG = HK*covKalman(:,:,k)*HK' + Rmeas;
        KG = covKalman(:,:,k)*HK'/(SG);
        xkalman(:,k) = xkalman(:,k) + KG*(zmes - HK*xkalman(:,k));
        covKalman(:,:,k) = (eye(dim) - KG*HK)*covKalman(:,:,k);               
    end
    kend = k;
end
toc(time1)
%% Save Results
saveResults =1;
if (exist('saveResults','var'))
    caseNum = 1;
    while( exist([caseStr,'_',num2str(caseNum)],'dir'))
        caseNum = caseNum+1;
    end
    folderName = [caseStr,'_',num2str(caseNum)];
    mkdir(folderName)
    try
        save([folderName,'/case: ', caseStr, datestr(now, 'dd-mmm-yyyy'),'.mat'],'caseStr','x' ...
             ,'cov','covKalman','dt','fFPE','op','x0','x0hat','xkalman','xGT','bdim','t','q','expec','n'  ...
             ,'dim','gridT','n','bcon','','pk')
    catch
        save([folderName,'/case: ', caseStr, datestr(now, 'dd-mmm-yyyy'),'.mat'])
    end
    
end
%% Check Plots

% State Vector
afigure
plot(t,xGT,t,expec,'.')
xlabel('time(s)')
zlabel('position')
for i=1:dim
    legString{i}= strcat('Kalman ' , num2str(i)); 
end
for i=1:dim
    legString{dim+i}=strcat('FPE ' , num2str(i)); 
end
legend(legString)
for k=1:length(t)
   trCov(k) = det(cov(:,:,k)); 
   trCovKalman(k) = det(covKalman(:,:,k)); 
end

% Covariance
afigure
plot(t,trCov,t,trCovKalman,'.')
xlabel('time')
ylabel('cov det')
legend('FPE','Kalman')
if (exist('saveResults','var'))
   saveas(gcf,[folderName,'/covVsKalman.eps'])
   saveas(gcf,[folderName,'/covVsKalman.png'])
   %saveas(gcf,[folderName,'covVsKalman.emf'])
   saveas(gcf,[folderName,'/covVsKalman.fig'])
end

% Error
afigure
plot(t,sqrt(sum((expec-xGT).^2)),t,sqrt(sum((xkalman-xGT).^2)))
xlabel('time(s)')
ylabel('mean error')
title('State Norm Error for the double integrator with measurements')
legend('Tensor','Kalman')
if (exist('saveResults','var'))
   saveas(gcf,[folderName,'/errorVsKalman.png'])
   figure
   plot2DslicesAroundPoint( pk{1}, x0hat, gridT)
   saveas(gcf,[folderName,'/initialDistribution.png'])
   saveas(gcf,[folderName,'/initialDistribution.eps'])
   saveas(gcf,[folderName,'/initialDistribution.fig'])
end


%% animated figure
hf = figure;
save_to_file = 1;
if (save_to_file && exist('saveResults','var') )
    FPS = 15;  
    pause(1)
    str_title = ['Probability Density Function Evolution'];
    writerObj = VideoWriter([folderName,'/pdf_gaussian_.avi']);
    writerObj.FrameRate = FPS;
    myVideo.Quality = 100;
    set(hf,'Visible','on');
    open(writerObj);
    set(gcf,'Renderer','OpenGL'); %to save to file
    pause(2)
end
if dim==1
    hf = figure;
    hold on
    ht_pdf = pcolor(gridT{1},gridT{2},double(pk{1})');    
    ht_kf = scatter(xkalman(1,1) ,xkalman(2,1) );
    h = colorbar;
    %set(h, 'ylim', [0 0.25])
    %caxis([0 15])
    set(ht_pdf, 'EdgeColor', 'none');
    xlabel('x')
    ylabel('y')
    grid on
    axis ([bdim(1,:),bdim(2,:)])
    grid on

    
    for k = 2:10:kend
         set ( ht_pdf, 'CData', double(pk{k})' );
         set ( ht_kf, 'XData', xkalman(1,k),'YData', xkalman(2,k) );
         drawnow
         pause(1.0/1000);
         if (save_to_file )
             M = getframe(gcf); %hardcopy(hf, '-dzbuffer', '-r0');
             writeVideo(writerObj, M);
         end
    end

else
    hold on
    handleSlices = plot2DslicesAroundPoint( pk{1}, x0, gridT);
    %ht_kf = scatter3(handleSlices{1,2}, xkalman(1,1) ,xkalman(2,1),5 );
    for k = 2:1:length(t)
        plot2DslicesAroundPoint( pk{k}, expec(:,k), gridT, handleSlices)
        %set ( ht_kf, 'XData', xkalman(1,k),'YData', xkalman(2,k) );
        pause(1.0/1000);
        if (save_to_file )
            M = getframe(gcf); %hardcopy(hf, '-dzbuffer', '-r0');
            writeVideo(writerObj, M);
        end        
    end
end

if (save_to_file  )
    close(writerObj);
end
%%


%%
for k=1:kend
    if (~isempty(zkcompressed{k}))
        plot2DslicesAroundPoint( zkcompressed{k}, x0, gridT);
    end
end

%%
nLambda = zeros(length(t),1);
normLambda = zeros(length(t),1);
qLambda = zeros(length(t),1);

for k = 1:kend 
    nLambda(k) = length(pk{k}.lambda);
    normLambda(k) = sum(pk{k}.lambda);
    qLambda(k) = pk{k}.lambda(1)/normLambda(k);
end
afigure
plot(t,qLambda)
title('Ration between the first two \lambda coefficients')
xlabel('time(s)')
ylabel('\lambda_1/\lambda_2')


afigure
plotLambda = zeros(length(t),max(nLambda));
for k = 1:kend
    plotLambda(k,1:nLambda(k)) = pk{k}.lambda;
end
plot(t,plotLambda,'d')
title('\lambda Evolution with time')
xlabel('time(s)')
ylabel('\lambda')

