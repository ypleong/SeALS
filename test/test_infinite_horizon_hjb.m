% Test finite horizon hjb, backward euler
% July 14, 2017

% Note:

addpath(genpath('../src'));
clear all

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

%% Define dynamics

% Linear system
% For d = 1, it's always converges correctly.
% For d = 2, it's mostly found the correct eigenvalue. Domain needs to be
% big enough.
%
% d = 2;
% x = sym('x',[d,1]); %do not change
% n = [51;53; ];
% bdim = [-5 5;-5 5;];
% bcon = {{'d',0,0},{'d',0,0}};%{{'d',0,0},{'d',0,0},{'d',0,0},{'d',0,0}};
% als_options = {100,20,'average',1e-7,1e-12,0.01,15};
% als_variant = {10,20};
% tol_err_op = 1e-6;
% c_err = 1e-6;
% ninputs = 2;
% AA = randn(d,d);
% BB = eye(d,ninputs);
% QQ = diag(ones(d,1));
% RR = 1/2*diag(ones(ninputs,1));
% PP = care(AA,BB,QQ,RR)
% B = BB;
% G = B;
% f = AA*x;
% noise_cov = diag(ones(ninputs,1));
% q = x'*QQ*x;
% R = RR;
% lambda = 1;%noise_cov*R;


% Smooth 2D
% f = 2*[x(1)^5-x(1)^3-x(1)+x(1)*x(2)^4; x(2)^5-x(2)^3-x(2)+x(2)*x(1)^4];
% G = [x(1) 0; 0 x(2)];
% B = BB;
% noise_cov = diag([1 1]);
% q = x'*QQ*x;
% R = diag([1 1]);
% lambda = 1;%noise_cov*R;

% Simple Pendulum 1D  
% Note: Working example, can compare controlled and not controlled, 
% Converged c =  1.2320 (Correct)
% d = 1;
% x = sym('x',[d,1]); %do not change
% n = [100;];
% bdim = [-pi pi-(2*pi/n(1));];
% bcon = {{'p'}};
% bsca = [1 1;];
% als_options = {100,20,'average',1e-7,1e-12,0.01,15};
% als_variant = {10,20};
% tol_err_op = 1e-5;
% f = [sin(x(1));];
% G = [1];
% B = G;
% noise_cov = 1;
% q = x(1)^2;
% R = 1;
% lambda = noise_cov*R;

% Simple Pendulum  
% Note: Working example, can compare controlled and not controlled, 
% Converged c =  1.3406
% d = 2;
% x = sym('x',[d,1]); %do not change
% n = [100; 103;];
% bdim = [-pi pi-(2*pi/n(1));-5 5;];
% bcon = {{'p'},{'d',0,0}};
% bsca = [1 1; 1 1;];
% als_options = {100,20,'average',1e-7,1e-12,0.01,15};
% als_variant = {10,20};
% tol_err_op = 1e-5;
% f = [x(2); sin(x(1));];
% G = [0; 1];
% B = G;
% noise_cov = 1;
% q = x(1)^2+x(2)^2;
% R = 1;
% lambda = noise_cov*R;

% Inverted Pendulum
% Note: The solution is too sharp to work
% d = 2;
% x = sym('x',[d,1]); %do not change
% n = [200; 203;];
% bdim = [-3*pi 3*pi-(6*pi/n(1));-10 10;];
% bcon = {{'p'},{'d',0,0}};
% bsca = [1 1; 1 1;];
% als_options = {100,20,'average',1e-7,1e-12,0.01,15};
% als_variant = {10,20};
% tol_err_op = 1e-5;
% m = 2; M = 8; l = .5; g = 9.8; mr = m/(m+M); den = 4/3-mr*cos(x(1))^2;
% f1 = (g/l)*sin(x(1))/den;
% f2 = -0.5*mr*x(2)^2*sin(2*x(1))/den;
% f = [x(2) ; f1+f2];
% G = [0.01 ; -mr/(m*l)*cos(x(1))/den];
% B = G;
% noise_cov = 10*pi;
% q = 0.1*x(1)^2+0.05*x(2)^2;
% R = 0.02;
% lambda = noise_cov*R;
% 
% % linearized dynamics
% AA = [0 1; (g/l)/(4/3-mr) 0];%randn(d,d);% 
% BB = [0.01; -mr/(m*l)/(4/3-mr)]; %eye(d);%
% QQ = diag([0.1 0.05]);
% RR = 1/2*R;
% PP = care(AA,BB,QQ,RR)

% Dubin's car
% d = 3;
% x = sym('x',[d,1]); %do not change
% n = [151; 153; 155];
% bdim = [-8 8; -8 8; -pi pi-(2*pi/(n(3)+1));];
% bcon = {{'d',0,0},{'d',0,0},{'p'}};
% bsca = ones(d,2);
% als_options = {100,20,'average',1e-7,1e-12,0.01,15};
% als_variant = {10,20};
% tol_err_op = 1e-5;
% c_err = 1e-8;
% 
% G = [0; 0; 1];
% B = G;
% uref = [0];
% f1 = [cos(x(3));sin(x(3)); 0];
% f = f1;
% noise_cov = diag(1);
% q = 10*(x(1).^2+x(2).^2);
% R = diag(1);
% lambda = 1;%noise_cov*R;

% linearized dynamics
% AA = [0 0 0; 0 0 1; 0 0 0 ];%randn(d,d);% 
% BB = [0; 0;1]; %eye(d);%
% QQ = diag(ones(d,1));
% RR = 1/2*R;
% PP = care(AA,BB,QQ,RR);

% VTOL
% Note: Works when no gravity, offset u by uref kind of work
d = 6;
x = sym('x',[d,1]); %do not change
n = [101; 103; 105; 107; 107; 101];
bdim = [-5 5; -5 5; -5 5; -5 5; -pi pi-(2*pi/(n(5)+1)); -5 5;];
bcon = {{'d',0,0},{'d',0,0},{'d',0,0},{'d',0,0},{'p'},{'d',0,0}};
bsca = ones(d,2);
als_options = {100,20,'average',1e-7,1e-12,0.01,15};
als_variant = {10,20};
tol_err_op = 1e-4;
c_err = 1e-6;

g = 9.8; eps = 0.01;
G = [0 0; -sin(x(5)) eps*cos(x(5)); 0 0; cos(x(5)) eps*sin(x(5)); 0 0; 0 1];
B = G;
uref = [g 0]';
f1 = [x(2); 0; x(4); -g; x(6); 0];
f2 = G*uref ;
f = f1 + f2;
noise_cov = diag([1 1]);
q = (1*x(3).^2 + 1*x(1).^2 + 1*x(5).^2+x(2).^2 +x(4).^2+x(6).^2);%x'*x;
R = diag([1 1]);
lambda = 1;%noise_cov*R;

% linearized dynamics
AA = [0 1 0 0 0 0; 0 0 0 0 -g 0; 0 0 0 1 0 0; zeros(1,d); 0 0 0 0 0 1; zeros(1,6);];%randn(d,d);% 
BB = [0 0; 0 eps; 0 0; 1 0; 0 0; 0 1]; %eye(d);%
QQ = diag(ones(d,1));
RR = 1/2*R;
PP = care(AA,BB,QQ,RR);

% Quadcopter
% d = 12;
% x = sym('x',[d,1]); %do not change
% n = 101*ones(d,1);
% bdim = [-3 3 ; -3 3; -3 3; -pi pi; -pi pi; -pi pi;...
%     -2 2; -2 2; -2 2; -pi pi-(2*pi/n(10)); -pi pi-(2*pi/n(11)); -pi pi-(2*pi/n(12))];
% bcon = { {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , ...
%     {'d',0,0} , {'d',0,0} , {'d',0,0} , {'p'} , {'p'} , {'p'} };
% bsca = ones(d,2);
% als_options = {100,20,'average',1e-7,1e-12,0.01,15};
% als_variant = {10,20};
% tol_err_op = 1e-6;
% 
% g = 9.8;
% hf1 = sin(x(12))*sin(x(10))+cos(x(12))*cos(x(10))*sin(x(11));
% hf2 = cos(x(12))*sin(x(11))*sin(x(10))-cos(x(10))*sin(x(12));
% hf3 = cos(x(11))*cos(x(12));
% G = [hf1 0 0 0; hf2 0 0 0; hf3 0 0 0; ...
%     0 1 0 0; 0 0 1 0; 0 0 0 1; ...
%     0 0 0 0; 0 0 0 0; 0 0 0 0; ...
%     0 0 0 0; 0 0 0 0; 0 0 0 0];
% B = G;
% uref = [g 0 0 0]';
% f2 = G*uref ;
% f1 = [0; 0; -g; 0; 0; 0; x(1); x(2); x(3); x(4); x(5); x(6)];
% f = f1+f2;
% ninputs = 4;
% noise_cov = eye(ninputs);
% q = x'*x;
% R = 0.01*eye(ninputs);
% lambda = 0.01;

% Burger PDE
% d = 12;
% x = sym('x',[d,1]); %do not change
% n = 51*ones(d,1);
% bsca = ones(d,2);
% bdim = [-1,1]; %[-2 2];
% bcon = {{'d', 0,0}}; % {{'v'}};
% region = [];
% for i = 2:d
%     bdim = [bdim;-1 1]; % [bdim; -2 2];
%     bcon{i} = {'d', 0, 0}; % {'v'};
% end
% als_options = {100,20,'average',1e-7,1e-12,0.01,15};
% als_variant = {10,20};
% tol_err_op = 1e-4;
% 
% gamma = 0.2;
% [A, AA, B] = spectralburgersdir(d, gamma, -0.5, -0.2);
% B(B == 0) = 0.01*ones(sum(B==0),1);
% G = B;
% f = A*x + x.* (AA*x);
% for ii = 1:d
%     f(ii) = expand(f(ii));
% end
% noise_cov = 6;
% q = x'*x;
% R = 4*gamma;
% lambda = noise_cov*R;


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

% for i=1:d
%     gridT{i} = linspace(bdim(i,1),bdim(i,2),n(i))';
%     nd(i) = n(i);
%     dxd(i) = abs(gridT{i}(2) - gridT{i}(1));
%     acc(:,i) = [2,2]';
% end
% 
% [D,D2,fd1,fd2] = makediffop(gridT,nd,dxd,acc,bcon,region);

[~,gridT,~,D,D2,fd1,fd2] = makediffopinfspectral(bdim,n,bcon,zeros(d,2));
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

%% Create boundary conditions

% [op] = makebcop(op,bcon,bsca,n,fd1);
% [op] = makebcopforward(op,bcon,n);
% [op] = makebcopspectral(op,bcon,bsca,n,fd1);

%% Create initial conditions
fprintf('Creating initial conditions ...\n');
start_PDEinit = tic;
U = cell(d,1);
for i = 1:d
    U{i} = %exp(-gridT{i}.^2*2);% matrandnorm(n(i),1);%
end
initTens = {ktensor(U)};

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
bc = initTens{1};
% [bc] = makebcbackward(initTens{1},bcon,n);

eigc(1) = 1/norm(initTens{1});
F_all{1} = initTens{1}*eigc(1);


iter_time = zeros(tend,1);
%%
for ind = 2:tend
    start_iter = tic;
    
%     [bc] = makebcforward(bcon,n);
%     F = SRMultV(op,F_all{ind-1});
%     F = F + norm(F)*bc;
    
%     [bc] = makebcbackward(F_all{ind-1},bcon,n) + makebcforward(bcon,n);
    bc = F_all{ind-1};
    [F, ~] = als_sys(op,bc,bc,tol_err_op,als_options,debugging, 0);
%     [F, ~] = als_sys_var(op,bc,bc,tol_err_op,als_options,als_variant,debugging, 0);
    
    if ncomponents(F) > ncomponents(F_all{ind-1})
        [F,~] = tenid(F,tol_err_op,1,9,'frob',[],fnorm(F),0);
        [F, ~] = als2(F,tol_err_op);
    end
    eigc(ind) = 1/norm(F);
    F_all{ind} = F*eigc(ind);
    iter_time(ind-1) = toc(start_iter);
    
    if mod(ind,10) == 0
        fprintf('Iteration: %d  Current tensor rank: %d \n', ind, ncomponents(F_all{ind}))
    end
    if abs(eigc(ind)-eigc(ind-1)) < c_err %|| norm(F_all{ind} - F_all{ind-1}) < tol_err_op
        fprintf('Iteration: %d  Current tensor rank: %d \n', ind, ncomponents(F_all{ind}))
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

figure;plot(eigc(1:ind))
xlabel('Iteration')
ylabel('Eigenvalue')
title(['Final value: ', num2str(eigc(ind)),' Correct linear value: ', num2str(trace(PP*BB*noise_cov*BB'))])
%%

figure
coord = zeros(d,1);
plot2DslicesAroundPoint( F_all{1}, coord, gridT,[],'surf');

%% Dimension = 1 

if d == 1
%% Plot solution
    figure;
    plot(gridT{1},double(F_all{ind}))
    xlabel('X(m)')
    ylabel('Desirability')
    
    
    
%% Plot result

    px = zeros(n,tend);
    for k=1:ind-1;%length(tt)
       px(:,k) = double(F_all{k}); 
    end
    figure
    hold on
    surf(gridT{1},1:ind-1,px(:,1:ind-1)','EdgeColor','none');
%     scatter3(ts_grid(:,1),ts_grid(:,2),ts_value,10,'filled');
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
%% Plot solution  
    figure
    surf(gridT{1},gridT{2},double(F_all{ind-1})', 'edgecolor','None');
    xlabel('x')
    ylabel('y')


end

%% Dimension > 2

if d > 2
    
%%    
%     Fori = EvalT(F_all{ind-1},zeros(d,1),gridT);
    dim_plot = [1 2]+0;
    figure
    coord = ceil(n/2);
    plot2Dslice(F_all{ind},dim_plot,coord,gridT);%*(1/Fori)
    axis([bdim(dim_plot(1),:),bdim(dim_plot(2),:)])
    
%% LQR 
    dim_plot = [2 3]+0;
    [gridx, gridy] = meshgrid(gridT{dim_plot(1)},gridT{dim_plot(2)});
    lqr_values = reshape(exp(-(sum(([gridx(:) gridy(:)]*PP(dim_plot,dim_plot)).*[gridx(:) gridy(:)],2))/lambda)...
        ,n(dim_plot(2)),n(dim_plot(1)));
    figure; 
    surf(gridT{dim_plot(1)},gridT{dim_plot(2)},lqr_values,'edgecolor','none');  
    xlabel('x1')
    ylabel('x2')
    axis([bdim(dim_plot(1),:),bdim(dim_plot(2),:)])
    
end


%% Simulations

saveplots = 0;
savedata = 0;

h = 0.0005;
[fFunc,GFunc,BFunc,noise_covFunc,qFunc,RFunc] = makefuncdyn(f1,G,B,noise_cov,q,R,x);
% sim_config = {20,h,repmat([1;1;-1;],1,1),[],[]};
sim_config = {20,h,repmat([0;0;1;0;1;0;],1,1),[],[]};
% sim_config = {20,h,repmat([1],1,1),[],[]};
sim_data = {lambda,gridT,R,noise_cov,F_all{ind},D,uref,fFunc,GFunc,BFunc,qFunc,bdim,bcon,region};

rng(100);
sim_run(sim_config,sim_data,saveplots,savedata,run,dirpath);


%% LQR Simulations

rng(100);
sim_data2 = {AA,BB,RR,QQ};
sim_run_lqg(sim_config,sim_data,sim_data2,saveplots,savedata,run,dirpath);
