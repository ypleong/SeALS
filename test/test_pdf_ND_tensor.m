
%% Test Tensor Form PDF Evolution ND
clear all


%% Create operator in tensor form
% Initialization

dim = 3;
x = sym('x',[dim,1]); %do not change
%n = [101,111,121,131,101,101];

% Initial distribution
sigma02 = 0.5;
for i=1:dim
    n(i) = 101+10*i;
    bdim(i,:) = [-5 25 ];
    x0(i) = 2+i/dim*10;
    diagSigma(i) = sigma02 + i/dim*0.3;
    dx(i) = (bdim(i,2)-bdim(i,1))/(n(i)-1);
    xvector{i} = [-5:dx(i):25]';
end


%p0 = zeros(n);
sigma02Matrix = diag(diagSigma); % [0.8 0.2; 0.2 sigma02];

for i=1:dim
   p0vector{i} =  normpdf(xvector{i},x0(i),sqrt(diagSigma(i)));
end


% szargs = cell( 1, dim  ); % We'll use this with ind2sub in the loop
% tic 
% for ii=1:numel(p0)
%     [ szargs{:} ] = ind2sub( size( p0 ), ii ); % Convert linear index back to subscripts
%     prodACC = 1;
%     for k=1:dim
%         prodACC = prodACC*p0vector{k}(szargs{k});
%     end
%     p0(szargs{:}) = prodACC;
% end
% toc

if dim == 2
    [xgrid,ygrid] = meshgrid(xvector,xvector);
    xyvector = [reshape(xgrid,n*n,1),reshape(ygrid,n*n,1)];
    p0 = reshape(mvnpdf(xyvector,x0,sigma02Matrix),n,n);
end
if dim == 2
    hh = pcolor(xvector{1},xvector{2},squeeze(p0(:,x0(2),:)));
    colorbar
    set(hh, 'EdgeColor', 'none');
    xlabel('x')
    zlabel('p(x)')
    grid on
end
    
% Create tensor structure
rankD = 10;

%[p0compressed] = cp_als(tensor(p0),rankD);


p0compressed = ktensor(p0vector);
pk{1} = p0compressed;

%% Tensor parameters

for i=1:dim
   bcon{i} = {'d',0,0};
end
bsca = []; %no manual scaling
region = [];
regval = 1;
regsca = [];
sca_ver = 1;

tol_err_op = 1e-6;
tol_err = 1e-9;
als_options = [];
als_variant = []; %{10,50};
debugging = 0;

run = 1;

fprintf(['Starting run ',num2str(run),' with main_run \n'])

for i=1:dim
    gridT{i} = xvector{i};
    nd(i) = n(i);
    dxd(i) = dx(i);
    acc(i,:) = [2,2]'
end
%bcon{1}{1} = 0;

[D,D2,fd1,fd2] = makediffop(gridT,nd,dxd,acc,bcon,region);


% Simulation Parameters
dt = 0.01;
finalt = 3;
t = 0:dt:finalt;
aspeed = zeros(dim,1);
qdiag = zeros(dim,1);
for i=1:dim
    qdiag(i) = 0.3 + i/dim*0.4;
    aspeed(i) = 4 - i/dim*3;
end
q =  diag( qdiag ); %[0.3,0.6]); %0.25;
%aspeed = [2;-2]; %[0,3]';

xkalman = zeros(dim,length(t));
covkalman = zeros(dim,dim,length(t));
expec = zeros(dim,length(t));
cov = zeros(dim,dim,length(t));

% Step up variables
cov(:,:,1) = sigma02Matrix;
expec(:,1) = x0;

% Kalman variables .. Equation xdot = a*x + noise
xkalman(:,1) = x0;
covKalman(:,:,1) = sigma02Matrix;



%% Create Cell Terms
aTens  = fcell2ftens( fsym2fcell(sym(aspeed) ,x), gridT);
dtTens = fcell2ftens( fsym2fcell(sym(dt)     ,x), gridT);
qTens  = fcell2ftens( fsym2fcell(sym(q)      ,x), gridT);
dtTen = DiagKTensor(dtTens{1}); 

% Create Operator
op = create_FP_op ( aTens, qTens, D, D2, dtTen, gridT, tol_err_op);

% Create weights for integration
for i=1:dim
    for j = 1:dim
        if i == j
            we{i,j} = xvector{j};
        else
            we{i,j} = ones(n(j),1);
        end         
    end
end
for i=1:dim
    weones{i} = ones(n(i),1);
end
weCov = cell(dim,dim,dim);
for i=1:dim
    for j = 1:dim
        weCov(i,j,:) = weones;
        if i == j
            weCov{i,j,i} = xvector{j}.*xvector{j};
        else
            weCov{i,j,i} = xvector{i};
            weCov{i,j,j} = xvector{j};
        end         
    end
end

%% Iterate over 
time1 = tic;
for k = 2:length(t)
    
    pk{k} = SRMultV( op, pk{k-1});
    [pk{k}, err_op, iter_op, enrich_op, t_step_op, cond_op, noreduce] = als2(pk{k},tol_err_op);
    
    for i=1:dim
        expec(i,k)  = intTens(pk{k}, [], gridT, we(i,:));
    end
    for i=1:dim
        for j = 1:dim
            cov(i,j,k) = intTens(pk{k}, [], gridT, weCov(i,j,:)) - expec(i,k)*expec(j,k);
        end
    end
    
    xkalman(:,k) = xkalman(:,k-1) + aspeed*dt;
    covKalman(:,:,k) = covKalman(:,:,k-1) + q*dt;
    
end
toc(time1)
%% Check Results

plot2d(xvector,double(pk{end}))
pcol
% Check Kalman 
figure
plot(t,xkalman,t,expec,'.')
xlabel('time')
zlabel('position')
%legString = [];
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
figure
plot(t,trCov,t,trCovKalman,'.')
xlabel('time')
zlabel('cov trace')
legend('FPE', 'Kalman')


%Compare 2D map kalman and Tensor-FKE
figure
plot2d( xvector, reshape(mvnpdf(xyvector, xkalman(:,k)',covKalman(:,:,end)),n,n)-double(pk{end}) )
title('Difference with kalman')



%% Measure
zmes = [12,6]';
Rmes = diag([1 1]); %0.5*[3,1;1,2]*[3,1;1,2];
pz = reshape(mvnpdf(xyvector,zmes',Rmes),n(2),n(1))';
rankD = 2;
[zkcompressed] = cp_als(tensor(pz),rankD);

pkpos = HadTensProd(pk{length(t)},zkcompressed);
pkpos = pkpos *(1/  intTens(pkpos, [], gridT, weones));

HK = eye(dim);

% Kalman Measurement
KalmanG = covKalman(:,:,end)*HK'*inv(HK*covKalman(:,:,end)*HK'+Rmes);
xposKalman = xkalman(:,end) + KalmanG*(zmes - HK*xkalman(:,end));
covKalmanPos = (eye(dim)-KalmanG*HK)*covKalman(:,:,end);


%% animated figure
figure
ht_pdf = pcolor(xvector,xvector,double(pk{k})');
h = colorbar;
%set(h, 'ylim', [0 0.25])
caxis([0 0.25])
set(ht_pdf, 'EdgeColor', 'none');
xlabel('x')
ylabel('y')
grid on
axis ([0,25,0,25])
grid on
for k = 2:10:length(t)
     set ( ht_pdf, 'CData', double(pk{k})' );
     drawnow
     pause(1.0/10.0);
end

%% animated figure difference
figure
ht_pdf = pcolor(xvector,xvector,reshape(mvnpdf(xyvector, xkalman(:,k)',covKalman(:,:,end)),n,n)-double(pk{end})');
h = colorbar;
%set(h, 'ylim', [0 0.25])
%caxis([0 0.25])
set(ht_pdf, 'EdgeColor', 'none');
xlabel('x')
ylabel('y')
grid on
axis ([0,25,0,25])
grid on
for k = 2:10:length(t)
     set ( ht_pdf, 'CData', reshape(mvnpdf(xyvector, xkalman(:,k)',covKalman(:,:,end)),n,n)-double(pk{k}));
     drawnow
     pause(1.0/10.0);
end
