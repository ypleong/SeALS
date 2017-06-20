% Test finite horizon hjb
% June 20, 2017

addpath(genpath('../src'));
clear all

%%
run = 1;
dirpath = ['./test_run/test_finite_horizon_hjb_euler/run_',num2str(run),'/'];
if 7~=exist(dirpath,'dir') 
    mkdir(dirpath); 
end
diary([dirpath,'run_',num2str(run),'output'])
fprintf('--------------------- New Run --------------------- \n\n')
disp(datetime)

start_whole = tic;

%% Initialization

d = 1;
x = sym('x',[d,1]); %do not change
n = 201; 

h = 0.005;%pi/(2*n); % time step size
tend = 8;
tt = 0:h:tend;
% tplot = 0.15;
% plotgap = round(tplot/h);
% h = tplot/plotgap;
% nplots = round(tend/tplot);

bdim = [-10 10];
bcon = { {'d', 0, 0} };
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

%% Define dynamics

lambda = 1;
% f = 2*[x(1)^5-x(1)^3-x(1)+x(1)*x(2)^4; x(2)^5-x(2)^3-x(2)+x(2)*x(1)^4]; % dynamics
% linearized dynamics
AA = -1; 
BB = 1; 
QQ = eye(d);
PP = care(AA,BB,QQ);

% full dynamics
f = AA*x; 
G = BB;
B = BB;

noise_cov = 1;
q = x'*QQ*x;
R = 1;

%% Calculate differentiation operators and finite difference matrices

fprintf('Creating differential operators ...\n');
start_diff = tic;

for i=1:d
    grid{i} = linspace(bdim(i,1),bdim(i,2),n)';
    nd(i) = n;
    dxd(i) = abs(grid{i}(2) - grid{i}(1));
    acc(:,i) = [2,2]';
end

[D,D2,fd1,fd2] = makediffop(grid,nd,dxd,acc,bcon,region);

% [n,grid,~,D,D2,fd1,fd2] = makediffopspectral(bdim,n,bcon,[]);
toc(start_diff)
end_diff = toc(start_diff);

%% Calculate dynamics
fprintf('Creating dynamics operator ...\n');
start_dynamics = tic;
[fTens,GTens,BTens,noise_covTens,qTens,RTens] = maketensdyn(f,G,B,noise_cov,q,R,x,grid);
toc(start_dynamics)
end_dynamics = toc(start_dynamics);

%% Calculate operator
fprintf('Creating PDE operator ...\n');
start_PDEop = tic;
[op,conv,diff] = makeop(fTens,BTens,noise_covTens,qTens,D,D2,0,lambda);
op = op*h + DiagKTensor(oneTens(d, n));
toc(start_PDEop)
end_PDEop = toc(start_PDEop);

%% Create boundary conditions

% create scaling for bc
if isempty(bsca) == 1 || isempty(regsca) == 1
    
    if sca_ver == 1
        [bscat, ~] = make_bc_sca_var(op,grid,region,bcon);    
    elseif sca_ver == 2
        [bscat, ~] = make_bc_sca(op,bcon,region,regval,als_options,fd1,grid,x,n);
    else
        error('wrong specification on boundary scaling');
    end
end

if isempty(bsca)
    bsca = bscat;
end

% make bc
% [op] = makebcop(op,bcon,bsca,n,fd1);
% [op] = makebcopspectral(op,bcon,bsca,n,fd1);

%% Create initial conditions
fprintf('Creating initial conditions ...\n');
start_PDEinit = tic;
init = exp(-x'*PP*x/lambda);
initTens  = fcell2ftens( fsym2fcell(sym(init) ,x), grid);
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

F_all = cell(1,tend+1);
F_all{1} = initTens{1};

tdata = 0;
data = double(initTens{1})';
t = 0;
% data_true = exp(-100*atan(tan(sqrt(6)*t/5+atan(sqrt(6)*tan(1-grid{1})))/sqrt(6)).^2)';
% iter_time = zeros(1, nplots*plotgap);

for ind = 2:length(tt)
%     fprintf([num2str(i) ': '])
%     for iter = 1:plotgap
        start_iter = tic;
%         fprintf([num2str(iter) ' '])
%         ind = ((i-1)*plotgap+iter+1);
%         t = t + h;
%         if isempty(als_variant) %original    
%             [F, err, iter, Fcond, enrich, t_step, illcondmat, maxit, maxrank, F_cell, B_cell, b_cell] = ...
%                 als_sys(op,bc,F_all{ind-1},tol_err,als_options,debugging, 0);
%             restart = []; %no restarts for original
%         else %variant
%             [F, err, iter, Fcond, enrich, t_step, illcondmat, maxit, maxrank, F_cell, B_cell, b_cell, restart] = ...
%                 als_sys_var(op,bc,F_all{ind-1},tol_err,als_options,als_variant,debugging, 0);
%         end
        F = SRMultV(op,F_all{ind-1});
        [F_all{ind}, ~] = als2(F,tol_err_op);
        
%         bc = F2;
%         F_all{ind} = F2; 
        iter_time(ind-1) = toc(start_iter);
%     end
%     fprintf('\n')
%     tdata = [tdata; t];
%     data(i+1,:) = double(F_all{ind})';
%     data_true(i+1,:) = exp(-100*atan(tan(sqrt(6)*t/5+atan(sqrt(6)*tan(1-grid{1})))/sqrt(6)).^2);
end

time_solve = toc(start_solve);
toc;
disp('Solution complete');

time_whole = toc(start_whole);

fprintf(['Run ',num2str(run),' with main_run is complete \n'])

diary off

%% Visualize results

% arrange input data
% F = arrange(F);
% plotsolve = {F,err,enrich,t_step,Fcond,grid};
% plotcomp = {op,err_op,enrich_op,t_step_op,cond_op};
% plotdebug = {F_cell,b_cell,B_cell};

%% Dimension = 1 

if d == 1
%% True solution

    xx_nn = 47;
    ts_grid = csvread([dirpath,'grid.csv']);
    x_grid = ts_grid(1:xx_nn:end,1);
    t_grid = ts_grid(1:xx_nn,2);
    ts_value = csvread([dirpath,'value.csv']);
    ts_value_square = reshape(ts_value,xx_nn,size(ts_grid,1)/xx_nn)';

%% Plot result

    px = zeros(n,length(tt));
    for k=1:length(tt)
       px(:,k) = double(F_all{k}); 
    end

    plottend = length(tt);
    figure
    hold on
    surf(grid{1},tt(end:-1:end-plottend+1),px(:,1:plottend)','EdgeColor','none');
    scatter3(ts_grid(:,1),ts_grid(:,2),ts_value,10,'filled');
    zlim([-0.2 1.2])
    % ylim([7.8 8])
    xlabel('X(m)')
    ylabel('Time(s)')
    zlabel('Desirability(x,t)')
    title('Desirability function evolution')
    colorbar
end


