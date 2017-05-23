
%% Test Tensor Form PDF Evolution
clear all



%% Create operator in tensor form
% Initialization

d = 1;
x = sym('x',[d,1]); %do not change
n = 201;

% Initial distribution
bdim = [-5 25];
dx = (bdim(2)-bdim(1))/(n-1);
x0 = 2;
sigma02 = 0.5;
xvector = [-5:dx:25]';
p0 = normpdf(xvector,x0,sqrt(sigma02));

% Tensor parameters
bcon = { {'d',0,0} };
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

[n,gridT,region,D,D2,fd1,fd2] = makediffopspectral(bdim,n,bcon,region);
gridT{1} = xvector;
bcon{1}{1} = 0;
[D,D2,fd1,fd2] = makediffop(gridT,n,dx,[2,2]',bcon,region);


% Simulation Parameters
dt = 0.01;
finalt = 15;
t = 0:dt:finalt;
px = zeros(length(xvector),length(t));
xkalman = zeros(length(t),1);
covkalman = zeros(length(t),1);
exp = zeros(length(t),1);
exp2 = zeros(length(t),1);
cov = zeros(length(t),1);
q = 0.25; %0.25;
aspeed = 1;




q = 0.25;
acell = {{{@(x1)aspeed, [1]}}};
qcell = {{{@(x1)q/2, [1]}}};
dtcell = {{{@(x1)dt, [1]}}};
aTens = fcell2ftens(acell,gridT);
qTens = fcell2ftens(qcell,gridT);
dtTens = fcell2ftens(dtcell,gridT);

icell = {{{@(x1)ones(size(x1)),1}}};
iTens = fcell2ftens(icell,gridT);
qTen = DiagKTensor(qTens{1}); 
aTen = DiagKTensor(aTens{1}); 
iTen = DiagKTensor(iTens{1}); 
dtTen = DiagKTensor(dtTens{1}); 
aop = -SRMultM(aTen,D{1}) + SRMultM(qTen,D2{1});
op = iTen + SRMultM(dtTen, aop); %DiagKTensor(kTens{1}) + D2{1} + D2{2};

%% Compress operator
op_uncomp = op; %save uncompressed op

fprintf('Attempt to compress operator, rank(op)=%d\n', ncomponents(op));
rank_op_uncomp = ncomponents(op);
tic;
start_compress = tic;

[op, err_op, iter_op, enrich_op, t_step_op, cond_op, noreduce] = als2(op,tol_err_op);

compress_time = toc(start_compress);
toc;
rank_op_comp = ncomponents(op);
fprintf('Number of components after compression, %d\n', ncomponents(op));

%debug operator
D2matrix = reshape(double(aop),n,n);
h_op = pcolor(D2matrix(20:30,20:30))
set(h_op,'edgecolor','none')
colorbar

% Kalman variables
xkalman(1) = x0;
covKalman(1) = sigma02;
exp(1) = x0;
cov(1) = sigma02;
px(:,1) = p0;

% Create tensor structure
p0cell{1} = p0;

% Iterate over 
p = ktensor( p0cell );
pk{1} = p;

for i = 2:length(t)
    
    pk{i} = SRMultV( op, pk{i-1});
    ss = sum(pk{i}.U{1});
    pk{i}.lambda = pk{i}.lambda*ss;
    pk{i}.U{1} = pk{i}.U{1}/ss;
    px(:,i) = double(pk{i});
    %get expected values
    exp(i)  = px(:,i)'*xvector*dx;
    exp2(i) = px(:,i)'*(xvector.^2)*dx;
    cov(i) = exp2(i) - exp(i).^2;
    cov(i) = sum(px(:,i).*(xvector - exp(i)).*(xvector - exp(i)))*dx;
    xkalman(i) = xkalman(i-1) + aspeed*dt;
    xkalman(i);
    covKalman(i) = covKalman(i-1) + q*dt;
    px(:,i) = px(:,i)/sum(px(:,i))/dx;
    pk{i} = pk{i}*(1/sum(px(:,i))/dx);
end

figure
plot(xvector,px(:,end),xvector,normpdf(xvector,xkalman(end),sqrt(covKalman(end))))
legend('px','Kalman')

% Measure
zmes = 11;
Rmes = 2;
zkcell{1} =   normpdf(xvector,zmes,sqrt(Rmes));
zk = ktensor( zkcell );

pkpos = HadTensProd(pk{length(t)},zk);

%pkpos = HadTensProd(tensor(ones(n,1)),pkpos)

% Plotting
plot(xvector, double(pk{length(t)}), xvector, double(zk), xvector, double(pkpos)/sum(double(pkpos))/dx)
legend('prior','measurement','posterior');

KalmanG = covKalman(end)/(covKalman(end)+Rmes);
xposKalman = xkalman(end) + KalmanG*(zmes - xkalman(end));
covKalmanPos = covKalman(end)*(1-KalmanG);


xpos = (double(pkpos)/sum(double(pkpos))/dx)'*xvector*dx;
covpos = sum(double(pkpos)/sum(double(pkpos))/dx.*(xvector - xpos).*(xvector - xpos))*dx;

% Print results
fprintf(['Measure update: Mean values ',num2str(xposKalman,4),' Kalman, ',num2str(xpos, 4) , ' Tensor-FPE \n'])
fprintf(['Measure update: Mean cov ',num2str(covKalmanPos,4),' Kalman, ',num2str(covpos, 4) , ' Tensor-FPE \n'])

%% Plot time evolution
figure
h=pcolor(xvector,t,px') %,'EdgeColor','none')
set(h, 'EdgeColor', 'none');
xlabel('X(m)')
ylabel('Time(s)')
zlabel('Pdf(x,t)')
title('Pdf Evolution')
colorbar

