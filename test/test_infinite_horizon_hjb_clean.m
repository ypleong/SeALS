% Test infinite horizon hjb
% Sept 24, 2017

% Note:

clear all
addpath(genpath('/home/yleong/tensor_toolbox'));
addpath(genpath('../src'));

%%
run = 1;
dirpath = ['./test_run/test_infinite_horizon_hjb/run_',num2str(run),'/'];
if 7~=exist(dirpath,'dir') 
    mkdir(dirpath); 
end
diary([dirpath,'run_',num2str(run),'output'])
fprintf('--------------------- New Run --------------------- \n\n')
disp(datetime)

start_whole = tic;

dynamic_choice = 'SimplePendulum';
discretization = 'spectral'; % 'uniform'; %
initialization = 'linear'; % 'simple', 'random'
saveF = 0;

%% Define dynamics

% dynamic choices

if strcmp(dynamic_choice,'VTOL')
    
    % c = 4.5619
    % c = 4.4794

    d = 6;
    ninputs = 2;
    x = sym('x',[d,1]); %do not change
    
    n = [101; 103; 105; 107; 107; 101]; %
    bdim = [-5 5; -5 5; -5 5; -5 5; -pi pi-(2*pi/(n(5)+1)); -5 5;];
    bcon = {{'d',0,0},{'d',0,0},{'d',0,0},{'d',0,0},{'p'},{'d',0,0}};
    bsca = ones(d,2);
    als_options = {100,25,'average',1e-7,1e-12,0.01,15};
    als_variant = {10,20};
    tol_err_op = 1e-6;
    c_err = 1e-8;
        
    g = 9.8; eps = 0.01;
    G = [0 0; -sin(x(5)) eps*cos(x(5)); 0 0; cos(x(5)) eps*sin(x(5)); 0 0; 0 1];
    B = G;
    uref = [g 0]';
    f1 = [x(2); 0; x(4); -g; x(6); 0];
    f2 = G*uref ;
    f = f1 + f2;
    noise_cov = diag([1 1]);
    q = x'*x;
    R = diag([1 1]);
    lambda = 1;%noise_cov*R;

    % linearized dynamics
    AA = [0 1 0 0 0 0; 
          0 0 0 0 -g 0; 
          0 0 0 1 0 0; 
          0 0 0 0 0 0; 
          0 0 0 0 0 1; 
          zeros(1,6);]; 
    BB = [0 0; 0 eps; 0 0; 1 0; 0 0; 0 1]; 
    QQ = diag(ones(d,1));
    RR = 1/2*R;
    PP = care(AA,BB,QQ,RR);
    trace(PP*BB*noise_cov*BB')

elseif strcmp(dynamic_choice,'Quadcopter')
 
    d = 12;
    ninputs = 4;
    x = sym('x',[d,1]); %do not change
    n = 151*ones(d,1);
    bdim = [-5 5 ; -5 5; -5 5; -pi pi; -pi pi; -pi pi;...
        -5 5; -5 5; -5 5; -pi pi-(2*pi/n(10)); -pi pi-(2*pi/n(11)); -pi pi-(2*pi/n(12))];
    bcon = { {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , ...
        {'d',0,0} , {'d',0,0} , {'d',0,0} , {'p'} , {'p'} , {'p'} };
    bsca = ones(d,2);
    als_options = {100,20,'average',1e-7,1e-12,0.01,15};
    als_variant = {10,20};
    tol_err_op = 1e-5;
    c_err = 1e-25;

    g = 9.8;
    hf1 = sin(x(12))*sin(x(10))+cos(x(12))*cos(x(10))*sin(x(11));
    hf2 = cos(x(12))*sin(x(11))*sin(x(10))-cos(x(10))*sin(x(12));
    hf3 = cos(x(11))*cos(x(12));
    G = [hf1 0 0 0; hf2 0 0 0; hf3 0 0 0; ...
        0 1 0 0; 0 0 1 0; 0 0 0 1; ...
        0 0 0 0; 0 0 0 0; 0 0 0 0; ...
        0 0 0 0; 0 0 0 0; 0 0 0 0];
    B = G;
    uref = [g 0 0 0]';
    f2 = G*uref ;
    f1 = [0; 0; -g; 0; 0; 0; x(1); x(2); x(3); x(4); x(5); x(6)];
    f = f1+f2;
    noise_cov = eye(ninputs);
    q = x'*x;
    R = eye(ninputs);
    lambda = 1;

    AA = [0 0 0 0 0 0 0 0 0 0 g 0; 
          0 0 0 0 0 0 0 0 0 0 0 -g; 
          zeros(4,d); 
          1 0 0 0 0 0 0 0 0 0 0 0; 
          0 1 0 0 0 0 0 0 0 0 0 0;
          0 0 1 0 0 0 0 0 0 0 0 0;
          0 0 0 1 0 0 0 0 0 0 0 0;
          0 0 0 0 1 0 0 0 0 0 0 0;
          0 0 0 0 0 1 0 0 0 0 0 0];
    BB = [zeros(2,ninputs); 1 0 0 0; 
          0 1 0 0; 0 0 1 0; 0 0 0 1;
          zeros(6,ninputs)]; 
    QQ = diag(ones(d,1));
    RR = 1/2*R;
    PP = care(AA,BB,QQ,RR);
    trace(PP*BB*noise_cov*BB')

elseif strcmp(dynamic_choice,'Linear')

    d = 4;
    ninputs = 3;
    x = sym('x',[d,1]); %do not change
    n = 51;
    bdim = [-5,5];
    bcon = {{'d', 0,0}};
    bsca = ones(d,2);
    for i = 2:d
        n = [n; n(i-1)+2];
        bdim = [bdim;-5 5]; 
        bcon{i} = {'d', 0, 0}; 
    end
    als_options = {100,20,'average',1e-7,1e-12,0.01,15};
    als_variant = {10,20};
    tol_err_op = 1e-5;
    c_err = 1e-7;

    AA = randn(d,d);
    load([dirpath,'run_linear4D.mat'],'AA');
    BB = eye(d,ninputs);
    QQ = diag(ones(d,1));
    RR = 1/2*diag(ones(ninputs,1));
    PP = care(AA,BB,QQ,RR);
    noise_cov = diag(ones(ninputs,1));
    trace(PP*BB*noise_cov*BB')
    
    B = BB;
    G = B;
    f = AA*x;
    f1 = f;
    uref = [0 0 0]';
    q = x'*QQ*x;
    R = 2*RR;
    lambda = 1;%noise_cov*R;
    
elseif strcmp(dynamic_choice,'SimplePendulum')
    
    d = 2;
    ninputs = 1;
    x = sym('x',[d,1]); %do not change
    n = [101; 103;];
    bdim = [-pi pi-(2*pi/(n(1)+1)); -5 5;];
    bcon = {{'p'},{'d',0,0}};
    bsca = ones(d,2);
    als_options = {100,25,'average',1e-7,1e-12,0.01,15};
    als_variant = {10,20};
    tol_err_op = 1e-6;
    c_err = 1e-7;

    % Simple Pendulum
    f = [x(2) ; sin(x(1))];
    f1 = f;
    uref = 0;
    G = [0.01 ; 1];
    B = G;
    noise_cov = 1;
    q = 0.1*x(1)^2+0.05*x(2)^2;
    R = 0.02;
    lambda = noise_cov*R;
    
    AA = [0 1; 
          1 0]; 
    BB = B;
    QQ = diag([0.1,0.05]);
    RR = 1/2*R;
    PP = care(AA,BB,QQ,RR);

elseif strcmp(dynamic_choice,'InvertedPendulum')
    
    d = 2;
    ninputs = 1;
    x = sym('x',[d,1]); %do not change
    n = [101; 103;];
    bdim = [-pi pi-(2*pi/(n(5)+1)); -10 10;];
    bcon = {{'p'},{'d',0,0}};
    bsca = ones(d,2);
    als_options = {100,25,'average',1e-7,1e-12,0.01,15};
    tol_err_op = 1e-6;
    c_err = 1e-7;
    
    m = 2; M = 8; l = .5; g = 9.8; mr = m/(m+M); den = 4/3-mr*cos(x(1))^2;
    ff1 = (g/l)*sin(x(1))/den;
    ff2 = -0.5*mr*x(2)^2*sin(2*x(1))/den;
    f = [x(2) ; ff1+ff2];
    f1 = f;
    G = [0.01 ; -mr/(m*l)*cos(x(1))/den];
    B = G;
    noise_cov = 10*pi;
    q = 0.1*x(1)^2+0.05*x(2)^2;
    R = 0.02;
    lambda = noise_cov*R;

    % linearized dynamics
    AA = [0 1; (g/l)/(4/3-mr) 0];%randn(d,d);% 
    BB = [0.01; -mr/(m*l)/(4/3-mr)]; %eye(d);%
    QQ = diag([0.1 0.05]);
    RR = 1/2*R;
    PP = care(AA,BB,QQ,RR);

end

%% Initialization

tend = 1000;

region = [];
regval = 1;
regsca = [];
sca_ver = 1;

debugging = 0;

fprintf(['Starting run ',num2str(run),' with main_run \n'])

%% Calculate differentiation operators and finite difference matrices

fprintf('Creating differential operators ...\n');
start_diff = tic;

if strcmp(discretization,'uniform')
    for i=1:d
        gridT{i} = linspace(bdim(i,1),bdim(i,2),n(i))';
        nd(i) = n(i);
        dxd(i) = abs(gridT{i}(2) - gridT{i}(1));
        acc(:,i) = [2,2]';
    end
    [D,D2,fd1,fd2] = makediffop(gridT,nd,dxd,acc,bcon,region);
elseif strcmp(discretization,'spectral')
    [~,gridT,~,D,D2,fd1,fd2] = makediffopinfspectral(bdim,n,bcon,zeros(d,2));
end

toc(start_diff)
end_diff = toc(start_diff);

%% Calculate dynamics
fprintf('Creating dynamics operator ...\n');
start_dynamics = tic;
[fTens,GTens,BTens,noise_covTens,qTens,RTens] = maketensdyn(f,G,B,noise_cov,q,R,x,gridT);
toc(start_dynamics)
end_dynamics = toc(start_dynamics);

%% Calculate operator
fprintf('Creating PDE operator ...\n');
start_PDEop = tic;
[op,conv,diff] = makeop(fTens,BTens,noise_covTens,qTens,D,D2,0,lambda);
op = -op*lambda;
toc(start_PDEop)
end_PDEop = toc(start_PDEop);

%% Create initial conditions
fprintf('Creating initial conditions ...\n');
start_PDEinit = tic;
U = cell(d,1);
numrank = 3;

if strcmp(initialization,'linear')

    for i = 1:d
        for j = i+1:d
            dim_plot = [i j];
            [gridx, gridy] = meshgrid(gridT{dim_plot(1)},gridT{dim_plot(2)});
            lqr_values = reshape(exp(-(sum(([gridx(:) gridy(:)]*PP(dim_plot,dim_plot)).*[gridx(:) gridy(:)],2))/lambda)...
                ,n(dim_plot(2)),n(dim_plot(1)))';   
            [fftest,~,out] = cp_als(tensor(lqr_values),numrank);
            if ~out.ill_cond
                U{dim_plot(1)} = [U{dim_plot(1)} bsxfun(@times,sqrt(fftest.lambda)',fftest.U{1})];
                U{dim_plot(2)} = [U{dim_plot(2)} bsxfun(@times,sqrt(fftest.lambda)',fftest.U{2})];
            else
                U{dim_plot(1)} = [U{dim_plot(1)} repmat(exp(-gridT{dim_plot(1)}.^2*PP(dim_plot(1),dim_plot(1))),1,numrank)];%ones(n(dim_plot(1)),numrank)];
                U{dim_plot(2)} = [U{dim_plot(2)} repmat(exp(-gridT{dim_plot(2)}.^2*PP(dim_plot(2),dim_plot(2))),1,numrank)];%ones(n(dim_plot(2)),numrank)];
            end
        end
    end
    
elseif strcmp(initialization,'simple')
    for i = 1:d
        U{i} = exp(-gridT{i}.^2*2);
    end
    
elseif strcmp(initialization,'random')
    for i = 1:d
        U{i} = matrandnorm(n(i),1);
    end
end

F = ktensor(U);
fprintf('Attempt to compress initial tensor, rank(Finit) = %d\n', ncomponents(F));
[F,~] = tenid(F,tol_err_op,1,9,'frob',[],fnorm(F),0);
fprintf('Number of components after TENID compression, %d\n', ncomponents(F));
[F, ~] = als2(F,tol_err_op);
fprintf('Number of components after ALS compression, %d\n', ncomponents(F));

initTens = {F};

toc(start_PDEinit)
end_PDEinit = toc(start_PDEinit);

%% Compress operator

op_uncomp = op; %save uncompressed op

fprintf('Attempt to compress operator, rank(op)=%d\n', ncomponents(op));
rank_op_uncomp = ncomponents(op);

start_compress_id = tic;
fprintf('Target CTD: %d terms above tol\n', length(find(op.lambda>tol_err_op)));
fprintf('Running TENID with frobenius norm:\n')
[op,~] = tenid(op,tol_err_op,1,9,'frob',[],fnorm(op),0);
op = fixsigns(arrange(op));
compress_time_id = toc(start_compress_id);
fprintf('Number of components after TENID compression, %d\n', ncomponents(op));
toc(start_compress_id)

start_compress = tic;
fprintf('Running ALS:\n')
[op, err_op, iter_op, enrich_op, t_step_op, cond_op, noreduce] = als2(op,tol_err_op);
rank_op_comp = ncomponents(op);
fprintf('Number of components after ALS compression, %d\n', ncomponents(op));
compress_time = toc(start_compress);
toc(start_compress)

%% Solve system

disp('Beginning Solving');
tic;
start_solve = tic;

F_all = cell(1,tend);
eigc = zeros(1,tend);

[bc] = initTens{1};

eigc(1) = 1/norm(initTens{1});
F_last = initTens{1}*eigc(1);

if saveF
    F_all{1} = F_last;
end

nndata = sum(n);
iter_time = zeros(tend,1);
errF = zeros(tend,1);


%%
for ind = 2:tend
    start_iter = tic;
    
    bc = F_last;
    [F, ~] = als_sys_var(op,bc,bc,tol_err_op,als_options,als_variant,debugging, 0);
    if ncomponents(F) > ncomponents(F_last)
        [F,~] = tenid(F,tol_err_op,1,9,'frob',[],fnorm(F),0);
        [F, ~] = als2(F,tol_err_op);
    end
    eigc(ind) = 1/norm(F);
    F = F*(eigc(ind));
    
    if saveF
    	F_all{ind} = F;
    end
    
    iter_time(ind-1) = toc(start_iter);
    
    F_diff = F - F_last;
    [F_diff,~] = tenid(F_diff,tol_err_op,1,9,'frob',[],fnorm(F_diff),0);
    
    errF(ind-1) = norm(F_diff)/(nndata*(ncomponents(F_diff)));
    
    F_last = F;
    
    if mod(ind,10) == 0
        fprintf('Iteration: %d  Current tensor rank: %d  Current eigenvalue: %.5f  Difference: %.5e\n',...
            ind, ncomponents(F),eigc(ind),errF(ind-1))
    end
    if errF(ind-1) < c_err
        
        fprintf('Iteration: %d  Current tensor rank: %d  Current eigenvalue: %.5f  Difference: %.5e\n\n',...
            ind, ncomponents(F),eigc(ind),errF(ind-1))
        break
    end
    
    save([dirpath,'run_',num2str(run),'data'])
end

time_solve = toc(start_solve);
toc;
disp('Solution complete');

time_whole = toc(start_whole);

fprintf(['Run ',num2str(run),' with main_run is complete \n'])

diary off

save([dirpath,'run_',num2str(run),'data'])

%% Plotting

dim_plot = [1 2];

%% Eigenvalues
figure;
plot(eigc(1:ind),'linewidth',1.5)
xlabel('Iteration')
ylabel('Eigenvalue')
title(['Final Value: ', num2str(eigc(ind)),' Linear Value: ', ...
    num2str(trace(PP*BB*noise_cov*BB'))])
set(gca,'fontsize',18)

%% Errors
figure;
plot(errF(1:ind),'linewidth',1.5)
xlabel('Iteration')
ylabel('Error')
set(gca,'fontsize',18)
% title(['Final error: ', num2str(errF(ind))])

%% Slices
figure
coord = zeros(d,1);
plot2DslicesAroundPoint(F_last, coord, gridT,[],'surf');

%%
figure;
coord = ceil(n/2);
dim_plot = [5 6];
plot2Dslice(F_last*(1/F_last(coord)),dim_plot,coord,gridT,[],[],'surf');

%% LQR 
[gridx, gridy] = meshgrid(gridT{dim_plot(1)},gridT{dim_plot(2)});
lqr_values = reshape(exp(-(sum(([gridx(:) gridy(:)]*PP(dim_plot,dim_plot)).*[gridx(:) gridy(:)],2))/lambda)...
    ,n(dim_plot(2)),n(dim_plot(1)));
figure; 
surf(gridT{dim_plot(1)},gridT{dim_plot(2)},lqr_values,'edgecolor','none');  
xlabel(['x',num2str(dim_plot(1))])
ylabel(['x',num2str(dim_plot(2))])
axis([bdim(dim_plot(1),:),bdim(dim_plot(2),:)])

%% Simulations

saveplots = 0;
savedata = 1;

h = 0.0005;
[fFunc,GFunc,BFunc,noise_covFunc,qFunc,RFunc] = makefuncdyn(f1,G,B,noise_cov,q,R,x);
% sim_config = {10,h,repmat([0;0;1;0],1,10),[],[]};
sim_config = {10,h,randi([-100, 100],d,10)/100,[],[]};
sim_data = {lambda,gridT,R,noise_cov,F_last,D,uref,fFunc,GFunc,BFunc,qFunc,bdim,bcon,region};
sim_data2 = {AA,BB,RR,QQ};

%%
rng(100);
sim_run(sim_config,sim_data,saveplots,savedata,run,dirpath);

%%
rng(100); 
sim_run_lqg(sim_config,sim_data,sim_data2,saveplots,savedata,run,dirpath);

