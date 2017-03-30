function test_run_spectral_helmholtz(run)

% Helmholtz equation (7.6) in Spectral Methods in MATLAB by Trefethen
% Dxx u + Dyy u + k^2 u = f(x,y)  -1 < x,y < 1
% u = 0 at boundary
% k = 9
% f(x,y) = e^(-10(y-1)^2) e^(-10(x-1/2)^2)

start_whole = tic;

%% Initialization

d = 2;
x = sym('x',[d,1]); %do not change
n = 101;

bdim = [-1 1 ; -1 1];
bcon = { {'d',0,0} , {'d',0,0} };
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

%% Calculate differentiation operators and finite difference matrices

[n,grid,region,D,D2,fd1,fd2] = makediffopspectral(bdim,n,bcon,region);


%% Calculate dynamics

%   x = sym('x',[2,1]) %two dimensions
%   f = x(1)*x(2)+x(2)
%   becomes in cell array represenation
%   fcell = { { @(x2)x2, [2] }, { @(x1)x1, [1] ; @(x2)x2, [2]} }

fcell = {{{@(x1)exp(-10*(x1-1/2).^2), [1]; @(x2)exp(-10*(x2-1).^2), [2]}}};
k = 9;
kcell = {{{@(x1)k^2, [1]}}};
fTens = fcell2ftens(fcell,grid);
kTens = fcell2ftens(kcell,grid);

op = DiagKTensor(kTens{1}) + D2{1} + D2{2};

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
[bc] = makebc(bcon,bsca,grid,x,n)+fTens{1};
[op] = makebcop(op,bcon,bsca,n,fd1);

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

if isempty(als_variant) %original    
    [F, err, iter, Fcond, enrich, t_step, illcondmat, maxit, maxrank, F_cell, B_cell, b_cell] = ...
        als_sys(op,bc,[],tol_err,als_options,debugging);
    restart = []; %no restarts for original
else %variant
    [F, err, iter, Fcond, enrich, t_step, illcondmat, maxit, maxrank, F_cell, B_cell, b_cell, restart] = ...
        als_sys_var(op,bc,[],tol_err,als_options,als_variant,debugging);
end

time_solve = toc(start_solve);
toc;
disp('Solution complete');

%% Visualize results

% arrange input data
F = arrange(F);
plotsolve = {F,err,enrich,t_step,Fcond,grid};
plotcomp = {op,err_op,enrich_op,t_step_op,cond_op};
plotdebug = {F_cell,b_cell,B_cell};

% plot results from run als and compress operator
try
    fprintf('Plotting results \n')
    visres(plotsolve,plotcomp,plotdebug,n,debugging,0,0,restart,run)
    fprintf('Plotting complete \n')
catch
    fprintf('Could not visualize results \n')
end

time_whole = toc(start_whole);

fprintf(['Run ',num2str(run),' with main_run is complete \n'])
