% Test finite horizon hjb, backward euler
% June 20, 2017

% Note:
% - Spectral work with backward euler, but in 2D ALS can be ill-conditioned
% - Domain needs to be really large to work properly for both finite
% difference and spectral

addpath(genpath('../src'));
clear all

%%
run = 4;
dirpath = ['./test_run/test_finite_horizon_hjb_backward_euler/run_',num2str(run),'/'];
if 7~=exist(dirpath,'dir') 
    mkdir(dirpath); 
end
diary([dirpath,'run_',num2str(run),'output'])
fprintf('--------------------- New Run --------------------- \n\n')
disp(datetime)

start_whole = tic;

%% Initialization

d = 2;
x = sym('x',[d,1]); %do not change
n = [201; 401;]; 

h = 0.001;%pi/(2*n); % time step size
tend = 2;
tt = 0:h:tend;

bdim = [-pi pi; -10 10];
bcon = { {'p'},{'d', 0, 0}};
bsca = []; %no manual scaling
region = [];
regval = 1;
regsca = [];
sca_ver = 1;

tol_err_op = 1e-6;
tol_err = 1e-9;
als_options = {100,20,'average',1e-7,1e-12,0.01,15};
als_variant = []; %{10,50};
debugging = 0;

fprintf(['Starting run ',num2str(run),' with main_run \n'])

%% Define dynamics

% full dynamics
% f = x.^3 + 5*x.^2 + x; 
% f = AA*x;

% Smooth 2D
% f = 2*[x(1)^5-x(1)^3-x(1)+x(1)*x(2)^4; x(2)^5-x(2)^3-x(2)+x(2)*x(1)^4];
% G = [x(1) 0; 0 x(2)];
% B = BB;
% noise_cov = diag([1 1]);
% q = x'*QQ*x;
% R = diag([1 1]);
% lambda = 1;%noise_cov*R;


% Inverted Pendulum
m = 2; M = 8; l = .5; g = 9.8; mr = m/(m+M); den = 4/3-mr*cos(x(1))^2;
f1 = (g/l)*sin(x(1))/den;
f2 = -0.5*mr*x(2)^2*sin(2*x(1))/den;
f = [x(2) ; f1+f2];
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
PP = care(AA,BB,QQ)/1000;


%% Calculate differentiation operators and finite difference matrices

fprintf('Creating differential operators ...\n');
start_diff = tic;

for i=1:d
    gridT{i} = linspace(bdim(i,1),bdim(i,2),n(i))';
    nd(i) = n(i);
    dxd(i) = abs(gridT{i}(2) - gridT{i}(1));
    acc(:,i) = [2,2]';
end

[D,D2,fd1,fd2] = makediffop(gridT,nd,dxd,acc,bcon,region);

% [~,gridT,~,D,D2,fd1,fd2] = makediffopspectral(bdim,n,bcon,[]);
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
op = DiagKTensor(oneTens(d, n)) + op*h;
toc(start_PDEop)
end_PDEop = toc(start_PDEop);

%% Create boundary conditions

% create scaling for bc
if isempty(bsca) == 1 || isempty(regsca) == 1
    
    if sca_ver == 1
        [bscat, ~] = make_bc_sca_var(op,gridT,region,bcon);    
    elseif sca_ver == 2
        [bscat, ~] = make_bc_sca(op,bcon,region,regval,als_options,fd1,gridT,x,n);
    else
        error('wrong specification on boundary scaling');
    end
end

if isempty(bsca)
    bsca = bscat;
end

% make bc
[op] = makebcop(op,bcon,bsca,n,fd1);
% [op] = makebcopspectral(op,bcon,bsca,n,fd1);

%% Create initial conditions
fprintf('Creating initial conditions ...\n');
start_PDEinit = tic;
if d == 2
    init = exp(-x(1)'*PP(1,1)*x(1)/lambda)*exp(-x(2)'*PP(2,2)*x(2)/lambda);
elseif d == 1
    init = exp(-x(1)'*PP(1,1)*x(1)/lambda);
end
initTens  = fcell2ftens( fsym2fcell(sym(init) ,x), gridT);
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

iter_time = zeros(length(tt),1);
%%
for ind = 421:length(tt)
    start_iter = tic;
%     [F, ~] = als_sys(op,F_all{ind-1},F_all{ind-1},tol_err_op,als_options,debugging, 0);
    F = SRMultV(op,F_all{ind-1});
    
    if ncomponents(F) > ncomponents(F_all{ind-1})
        [F,~] = tenid(F,tol_err_op,1,9,'frob',[],fnorm(F),0);
        [F_all{ind}, ~] = als2(F,tol_err_op);
    else
        F_all{ind} = F;
    end
    iter_time(ind-1) = toc(start_iter);
    
    if mod(ind,10) == 0
        fprintf('Time: %.4fs  Current tensor rank: %d \n', tt(ind), ncomponents(F_all{ind}))
    end
    if norm(F)/prod(n) < tol_err_op/10
        break
    end
end

time_solve = toc(start_solve);
toc;
disp('Solution complete');

time_whole = toc(start_whole);

fprintf(['Run ',num2str(run),' with main_run is complete \n'])

diary off

save([dirpath,'run_',num2str(run),'data'])
    
%% Visualize results

%% Dimension = 1 

if d == 1
%% True solution

    ts_grid = csvread([dirpath,'grid.csv']);
    ts_value = csvread([dirpath,'value.csv']);

%% Plot result

    px = zeros(n,length(tt));
    for k=1:ind-1;%length(tt)
       px(:,k) = double(F_all{k}); 
    end

    plottend = length(tt)-length(tt)+ind-1;
    figure
    hold on
    surf(gridT{1},tt(end:-1:end-plottend+1),px(:,1:plottend)','EdgeColor','none');
    scatter3(ts_grid(:,1),ts_grid(:,2),ts_value,10,'filled');
    zlim([-0.2 1.2])
%     ylim([2.5 3])
    xlabel('X(m)')
    ylabel('Time(s)')
    zlabel('Desirability(x,t)')
    title('Desirability function evolution')
    colorbar
end

%% Dimension = 2

if d == 2
%% True solutions for linear system
    
% AA = eye(d);
% BB = eye(d);
% QQ = eye(d);
% PP = diag([6 6]);
% noise = diag([1 1]);
% lambda = 1; 

    ts_grid = csvread([dirpath,'grid 2.csv']);
    ts_value = csvread([dirpath,'value 2.csv']);

%% Plot result
    figure
    ht_pdf = surf(gridT{1},gridT{2},double(F_all{1})');
    h = colorbar;
    %set(h, 'ylim', [0 0.25])
%     caxis([0 1])
    set(ht_pdf, 'EdgeColor', 'none');
    xlabel('x')
    ylabel('y')
    title(['Time = ',num2str(tt(1))]);
    grid on
    axis ([bdim(1,:),bdim(2,:)])
    grid on
    for k = 2:10:ind%length(tt)
         set(ht_pdf, 'ZData', double(F_all{k})' );
         title(['Time = ',num2str(tt(k))]);
         drawnow
         pause(1.0/1000);
    end
    
%%  
    timenn = 115;
    k = 1;
    figure
    hold on
    surf(gridT{1},gridT{2},double(F_all{end})','EdgeColor', 'none');
    scatter3(ts_grid(1:timenn:end,1),ts_grid(1:timenn:end,2),ts_value(k-1+(1:timenn:end)),10,'filled');
    colorbar;
    xlabel('x')
    ylabel('y')
    xlim([-5 5])
    ylim([-5 5])

%%
    timenn = 115;
    figure
    ht_true = scatter3(ts_grid(1:timenn:end,1),ts_grid(1:timenn:end,2),ts_value(1:timenn:end),10,'filled');
    h = colorbar;
    %set(h, 'ylim', [0 0.25])
    caxis([0 1])
    xlabel('x')
    ylabel('y')
    title(['Time = ',num2str(tt(1))]);
    grid on
    axis ([bdim(1,:),bdim(2,:)])
    grid on
    for k = 2:1:timenn
         set(ht_true,'ZData', ts_value(k-1+(1:timenn:end)));
         title(['Time = ',num2str(ts_grid(k,3))]);
         drawnow
         pause(1.0/10);
    end

end
