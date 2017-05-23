function test_run_spectral_heat()

% Wave equation (3.13, pg 24) in Spectral Methods in MATLAB by Trefethen
% Dt u - k Dxx u = 0 -3 < th < 3, t>0
% u(x,0) = delta(x)

% From mathematica
% True solution = exp(-x^2/4kt)/sqrt(4kt)

start_whole = tic;

%% Initialization

d = 2;
x = sym('x',[d,1]); %do not change
n = 101; % better be odd

h = pi/(2*n); % time step size
tend = 8;
tplot = 0.15;
plotgap = round(tplot/h);
h = tplot/plotgap;
nplots = round(tend/tplot);

bdim = [-10 10];
bcon = {{'d',0,0}};
for i = 2:d
    bdim = [bdim; -10 10];
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
als_variant = {inf,inf}; %{10,50};
debugging = 0;

run = 1;

fprintf(['Starting run ',num2str(run),' with main_run \n'])

%% Calculate differentiation operators and finite difference matrices

[n,grid,region,D,D2,fd1,fd2] = makediffopspectral(bdim,n,bcon,region);


%% Calculate dynamics

%   x = sym('x',[2,1]) %two dimensions
%   f = x(1)*x(2)+x(2)
%   becomes in cell array represenation
%   fcell = { { @(x2)x2, [2] }, { @(x1)x1, [1] ; @(x2)x2, [2]} }

k = 0.1;
ccell = {{{@(x1)ones(size(x1))*-k, 1}}};
cTens = fcell2ftens(ccell,grid);
cTen = DiagKTensor(cTens{1});

ucell = {{{@(x1)(exp(-x1.^2/(4*k))/sqrt(4*k)), 1}}};
for i = 2:d
    ucell{1,1,i} = {@(x1)(exp(-x1.^2/(4*k))/sqrt(4*k)), i};
end
uTens = fcell2ftens(ucell,grid);
hcell = {{{@(x1)ones(size(x1))*h,1}}};
hTens = fcell2ftens(hcell,grid);
icell = {{{@(x1)ones(size(x1)),1}}};
iTens = fcell2ftens(icell,grid);

prodkD = SRMultM(cTen, D2{1});
for i = 2:d
    prodkD = prodkD + SRMultM(cTen, D2{i});
end

op = SRMultM(DiagKTensor(hTens{1}), prodkD) + DiagKTensor(iTens{1});

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
[bc] = uTens{1};
[op] = makebcopspectral(op,bcon,bsca,n,fd1);

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

%% Step 11: solve system

disp('Beginning Solving');
tic;
start_solve = tic;

F_all = cell(1,tend+1);
F_all{1} = bc;

data = double(bc)';
t = 1;
tdata = t;
data_true = data;
data_err = max(abs(data(:) - data_true(:)));
iter_time = zeros(1, nplots*plotgap);

for i = 1:nplots
    fprintf([num2str(i) ': '])
    for iter = 1:plotgap
        start_iter = tic;
        fprintf([num2str(iter) ' '])
        ind = ((i-1)*plotgap+iter+1);
        t = t + h;
        if isempty(als_variant) %original    
            [F, err, iter, Fcond, enrich, t_step, illcondmat, maxit, maxrank, F_cell, B_cell, b_cell] = ...
                als_sys(op,bc,F_all{ind-1},tol_err,als_options,debugging, 0);
            restart = []; %no restarts for original
        else %variant
            [F, err, iter, Fcond, enrich, t_step, illcondmat, maxit, maxrank, F_cell, B_cell, b_cell, restart] = ...
                als_sys_var(op,bc,F_all{ind-1},tol_err,als_options,als_variant,debugging, 0);
        end
        [F2, err_op, iter_op, ~] = als2(F,tol_err_op);
        
        bc = F2;
        F_all{ind} = F2; 
        iter_time(ind-1) = toc(start_iter);
    end
    fprintf('\n')
    tdata = [tdata; t];
    if d ==1
        data = double(bc)';
        data_true = (exp(-grid{1}.^2/(4*k))/sqrt(4*k*t));
    elseif d == 2
        data = double(bc)';
        ind = ceil(n/2);
        data = bc(ind,ind);
        data_true = (1/sqrt(4*k*t));
    end
    data_err(ii+1) = max(abs(data(:) - data_true(:)));
end

time_solve = toc(start_solve);
toc;
disp('Solution complete');

%% Visualize results

% arrange input data
% F = arrange(F);
% plotsolve = {F,err,enrich,t_step,Fcond,grid};
% plotcomp = {op,err_op,enrich_op,t_step_op,cond_op};
% plotdebug = {F_cell,b_cell,B_cell};

% plot results from run als and compress operator
try
    fprintf('Plotting results \n')
%     visres(plotsolve,plotcomp,plotdebug,n,debugging,0,0,restart,run)
    if d == 1
    figure;
    hold on
    h1 = waterfall(grid{1},tdata,data);
    set(h1,'FaceColor','none')
    h2 = waterfall(grid{1},tdata,data_true);
    axis([-10 10 1 tend 0 5])
    ylabel('t')
    zlabel('u')
    hold off
    
    figure; 
    waterfall(grid{1},tdata,abs(data-data_true))
    ylabel('t')
    zlabel('abs(error)')
    
    else
        figure;
        plot(tdata,data_err)
        xlabel('t')
        ylabel('Error')
    end
    
    figure
    plot(iter_time)
    xlabel('Iteration')
    ylabel('Computation time (s)')
    
    fprintf('Plotting complete \n')
catch
    fprintf('Could not visualize results \n')
end

time_whole = toc(start_whole);

fprintf(['Run ',num2str(run),' with main_run is complete \n'])
