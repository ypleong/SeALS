
%% Test Tensor Form PDF Evolution ND
clear all



%% Create operator in tensor form
% Initialization

dim = 2; %only dim=2 here
x = sym('x',[dim,1]); %do not change
n = 301;

% Initial distribution
bdim = [-5 25 ; -5 25];
dx = (bdim(1,2)-bdim(1,1))/(n-1);
x0 = [7,13];
sigma02 = 0.5;
sigma02Matrix = [0.8 0.2; 0.2 sigma02];
xvector = [-5:dx:25]';
[xgrid,ygrid] = meshgrid(xvector,xvector);
xyvector = [reshape(xgrid,n*n,1),reshape(ygrid,n*n,1)];
p0 = reshape(mvnpdf(xyvector,x0,sigma02Matrix),n,n);

% Tensor parameters
bcon = { {'d',0,0}, {'d',0,0}};
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

if dim == 2
    hh = pcolor(xvector,xvector,p0);
    colorbar
    set(hh, 'EdgeColor', 'none');
    xlabel('x')
    zlabel('p(x)')
    grid on
end

fprintf(['Starting run ',num2str(run),' with main_run \n'])

for i=1:dim
    gridT{i} = xvector;
    nd(i) = n;
    dxd(i) = dx;
    acc(i,:) = [2,2]';
end

[D,D2,fd1,fd2] = makediffop(gridT,nd,dxd,acc,bcon,region);


% Simulation Parameters
dt = 0.01;
finalt = 5;
t = 0:dt:finalt;
q =  diag([0.3,0.6]); %0.25;
aspeed = [2;-2]; %[0,3]';

xkalman = zeros(dim,length(t));
covkalman = zeros(dim,dim,length(t));
expec = zeros(dim,length(t));
cov = zeros(dim,dim,length(t));


cov(:,:,1) = sigma02Matrix;
expec(:,1) = x0;

% Kalman variables .. Equation xdot = a*x + noise
xkalman(:,1) = x0;
covKalman(:,:,1) = sigma02Matrix;



%% Create Cell Terms
%acell = {{{@(x1)aspeed(1), [1]},{@(x2)aspeed(2), [2]}}};
%aTens = fcell2ftens(acell,gridT);

%dtcell = {{{@(x1)dt, [1]}}};
%dtTens = fcell2ftens(dtcell,gridT);
aTens  = fcell2ftens( fsym2fcell(sym(aspeed) ,x), gridT);
dtTens = fcell2ftens( fsym2fcell(sym(dt)     ,x), gridT);
qTens  = fcell2ftens( fsym2fcell(sym(q)      ,x), gridT);
dtTen = DiagKTensor(dtTens{1}); 

% Create Operator
op = create_op ( aTens, qTens, D, D2, dtTen, gridT, tol_err_op);

% Create tensor structure
rankD = 5;
[p0compressed] = cp_als(tensor(p0),rankD);

pk{1} = p0compressed;

%% Iterate over 
tic
for k = 2:length(t)
    
     pk{k} = SRMultV( op, pk{k-1});
    [pk{k}, err_op, iter_op, enrich_op, t_step_op, cond_op, noreduce] = als2(pk{k},tol_err_op);
    
%     we{1} = ones(1,n)*dx;
%     we{2} = xvector'*dx;
%     
    for i=1:dim
       we{i} = ones(1,n)*dx;
       wex{i} = xvector'*dx;
    end
%     intTens(pk{k},we, [1,2] )
%     intTens(pk{k},we, [2,1] )
%     expec(:,k) = [ intTens( intTens(pk{k},we, 2 ) ,wex,1)  
%                    intTens( intTens(pk{k},we, 1 ) ,wex,2)  ]
               
    %intTens(pk{k},we, [1,2] );
    expec(:,k)  = [ trapz( xvector, xvector'.*trapz(xvector,double(pk{k})',1 ) )
                  trapz( xvector, xvector.*trapz(xvector,double(pk{k})',2 ) )];
              
 
    xkalman(:,k) = xkalman(:,k-1) + aspeed*dt;
    covKalman(:,:,k) = covKalman(:,:,k-1) + q*dt;
    
end
toc
%% Check Results
plot2d(xvector,double(pk{end}))

% Check Kalman 
figure
plot(t,xkalman,t,expec,'.')
xlabel('time')
zlabel('position')
legend('KalmaX','KalmanY','FPE X','FPE Y' )

%Compare 2D map kalman and Tensor-FKE
figure
plot2d( xvector, reshape(mvnpdf(xyvector, xkalman(:,k)',covKalman(:,:,end)),n,n)-double(pk{end}) )
title('Difference with kalman')



%% Measure
zmes = [12,17]';
Rmes = diag([1 1]); %0.5*[3,1;1,2]*[3,1;1,2];
pz = reshape(mvnpdf(xyvector,zmes',Rmes),n,n);
rankD = 2;
[zkcompressed] = cp_als(tensor(pz),rankD);

pkpos = HadTensProd(pk{length(t)},zkcompressed);
pkpos = pkpos *(1/ intTens(pkpos,we, [1,2] ));

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
