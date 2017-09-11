% Finite Horizon
% Sept 7, 2017

addpath(genpath('../src'));
clear all

%%
run = 1;
dirpath = ['./test_run/test_finite_horizon_hjb_scaling/run_',num2str(run),'/'];
if 7~=exist(dirpath,'dir') 
    mkdir(dirpath); 
end
diary([dirpath,'run_',num2str(run),'output'])
fprintf('--------------------- New Run --------------------- \n\n')
disp(datetime)

start_whole = tic;

dynamic_choice = 'Linear';
time_stepping_choice = 'backward'; % forward

%% Define dynamics

% dynamic choices

if strcmp(dynamic_choice,'VTOL')

    d = 6;
    ninputs = 2;
    x = sym('x',[d,1]); %do not change
    
    n = [101; 103; 105; 107; 107; 101];
    bdim = [-5 5; -5 5; -5 5; -5 5; -pi pi-(2*pi/(n(5)+1)); -5 5;];
    bcon = {{'d',0,0},{'d',0,0},{'d',0,0},{'d',0,0},{'p'},{'d',0,0}};
    bsca = ones(d,2);
    als_options = {100,25,'average',1e-7,1e-12,0.01,15};
    als_variant = {10,20};
    tol_err_op = 1e-6;
    c_err = 1e-7;
        
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
          zeros(1,d); 
          0 0 0 0 0 1; 
          zeros(1,6);]; 
    BB = [0 0; 0 eps; 0 0; 1 0; 0 0; 0 1]; 
    QQ = diag(ones(d,1));
    RR = 1/2*R;
    PP = care(AA,BB,QQ,RR);

elseif strcmp(dynamic_choice,'Quadcopter')
 
    d = 12;
    ninputs = 4;
    x = sym('x',[d,1]); %do not change
    n = 51*ones(d,1);
    bdim = [-5 5 ; -5 5; -5 5; -pi pi; -pi pi; -pi pi;...
        -5 5; -5 5; -5 5; -pi pi-(2*pi/n(10)); -pi pi-(2*pi/n(11)); -pi pi-(2*pi/n(12))];
    bcon = { {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , ...
        {'d',0,0} , {'d',0,0} , {'d',0,0} , {'p'} , {'p'} , {'p'} };
    bsca = ones(d,2);
    als_options = {100,20,'average',1e-7,1e-12,0.01,15};
    als_variant = {10,20};
    tol_err_op = 1e-4;
    c_err = 1e-4;

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

elseif strcmp(dynamic_choice,'Linear')

    d = 2;
    ninputs = 1;
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
    BB = eye(d,ninputs);
    QQ = diag(ones(d,1));
    RR = 1/2*diag(ones(ninputs,1));
    PP = care(AA,BB,QQ,RR);
    B = BB;
    G = B;
    f = AA*x;
    f1 = f;
    noise_cov = diag(ones(ninputs,1));
    q = x'*QQ*x;
    R = RR;
    lambda = 1;%noise_cov*R;
    
elseif strcmp(dynamic_choice,'SimplePendulum')
    
    d = 2;
    ninputs = 1;
    x = sym('x',[d,1]); %do not change
    n = [101; 103;];
    bdim = [-pi pi-(2*pi/(n(5)+1)); -5 5;];
    bcon = {{'p'},{'d',0,0}};
    bsca = ones(d,2);
    als_options = {100,25,'average',1e-7,1e-12,0.01,15};
    tol_err_op = 1e-6;
    c_err = 1e-7;

    % Simple Pendulum
    f = [x(2) ; sin(x(1))];
    f1 = f;
    G = [0.01 ; 1];
    B = G;
    noise_cov = 1;
    q = 0.1*x(1)^2+0.05*x(2)^2;
    R = 0.02;
    lambda = noise_cov*R;
    
    AA = [0 1; 
          1 0]; 
    BB = B;
    QQ = diag(ones(d,1));
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

h = 0.001;  % time step size
tend = 5;
tt = 0:h:tend;

region = [];
regval = 1;
regsca = [];
sca_ver = 1;

debugging = 0;

fprintf(['Starting run ',num2str(run),' with main_run \n'])


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

% [~,gridT,~,D,D2,fd1,fd2] = makediffopspectral(bdim,n,bcon,[0 0]);
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
if strcmp(time_stepping_choice, 'forward')
    op = DiagKTensor(oneTens(d, n)) + op*h; % forward
elseif strcmp(time_stepping_choice, 'backward')
    op = DiagKTensor(oneTens(d, n)) - op*h; % backward
else
    error([time_stepping_choice, ' time stepping scheme does not exist.'])
end
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

% NOTE: This is bad for backward euler. Not if the boundary condition
% scaling is 1 for Dirichlet BC. For periodic and neumann, need to adjust 
% the boundary of F_all as well.
% [op] = makebcop(op,bcon,bsca,n,fd1);
% [op] = makebcopspectral(op,bcon,bsca,n,fd1);

if strcmp(time_stepping_choice, 'forward')
    [op] = makebcopforward(op,bcon,n); % forward
elseif strcmp(time_stepping_choice, 'backward')
    [op] = makebcop(op,bcon,bsca,n,fd1); % backward
else
    error([time_stepping_choice, ' time stepping scheme does not exist.'])
end

%% Create initial conditions
fprintf('Creating initial conditions ...\n');
start_PDEinit = tic;

if strcmp(dynamic_choice,'Linear')
    if d == 2
        [gridx, gridy] = meshgrid(gridT{1},gridT{2});
        init = reshape(exp(-(sum(([gridx(:) gridy(:)]*PP).*[gridx(:) gridy(:)],2))/lambda),n(2),n(1));
        initTens  = {cp_als(tensor(init'),10)};
    else
        init = exp(-x(1)'*PP*x(1)/lambda);
        initTens  = fcell2ftens( fsym2fcell(sym(init) ,x), gridT);
    end
else
    U = cell(d,1);
    for i = 1:d
        U{i} = exp(-gridT{i}.^2);
    end
    initTens = {ktensor(U)};
end
    
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

nt = length(tt);

F_all = cell(1,nt+1);
F_all{1} = initTens{1};
F_scale = ones(1,nt+1);

t = 0;

iter_time = zeros(nt,1);
%%
for ind = 2:nt+1
    start_iter = tic;
    if strcmp(time_stepping_choice, 'forward')
        [bc] = makebcforward(bcon,n);
        F = SRMultV(op,F_all{ind-1}) + bc;
    elseif strcmp(time_stepping_choice, 'backward')
        [bc] = makebcbackward(F_all{ind-1},bcon,n);
        [F, ~] = als_sys(op,bc,bc,tol_err_op,als_options,debugging, 0);
    else
        error([time_stepping_choice, ' time stepping scheme does not exist.'])
    end
    
    if ncomponents(F) > ncomponents(F_all{ind-1})
        [F,~] = tenid(F,tol_err_op,1,9,'frob',[],fnorm(F),0);
        [F, ~] = als2(F,tol_err_op);
    end
    
    F_norm = norm(F)/sqrt(ncomponents(F)*sum(n));
    F_scale(ind) = F_scale(ind-1)*F_norm;
    F_all{ind} = F*(1/F_norm);
    iter_time(ind-1) = toc(start_iter);
    
    if mod(ind,20) == 0
        fprintf('Time: %.4fs  Current tensor rank: %d \n', tt(ind), ncomponents(F_all{ind}))
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

%% Plot result

    px = zeros(n,length(tt));
    for k=1:ind-1;%length(tt)
       px(:,k) = double(F_all{k}*F_scale(k)); 
    end

    plottend = length(tt)-length(tt)+ind-1;
    figure
    hold on
    surf(gridT{1},tt(end:-1:end-plottend+1),px(:,1:plottend)','EdgeColor','none');
    zlim([0 1])
    xlabel('X(m)')
    ylabel('Time(s)')
    zlabel('Desirability(x,t)')
    title('Desirability function evolution')
    colorbar
end

%% Dimension > 1

if d == 2
    
%%    
    dim_plot = [1 2];
    coord = ceil(n/2);
    for k = 1:10:ind-1
         plot2Dslice(F_all{k}*F_scale(k),dim_plot,coord,gridT);
         colorbar;
%          caxis([0 1])
         axis ([bdim(dim_plot(1),:),bdim(dim_plot(2),:)])
         title(['Time = ',num2str(tt(k))]);
         zlim([0 100000])
         pause(1.0/1000);
    end
    
end


%% Simulations

saveplots = 0;
savedata = 0;
ninputs = 1;
uref = zeros(ninputs,1);

[fFunc,GFunc,BFunc,noise_covFunc,qFunc,RFunc] = makefuncdyn(f1,G,B,noise_cov,q,R,x);
sim_config = {h*(ind-1),h,repmat([0.2; 0; ],1,1),[],[]};
sim_data = {lambda,gridT,R,noise_cov,F_all,uref,D,fFunc,GFunc,BFunc,qFunc,bdim,bcon,region};

sim_finite_run(sim_config,sim_data,saveplots,savedata,run,dirpath)

%% LQR

SS = zeros(d,d,length(tt));
qq = zeros(length(tt),1);
SS(:,:,1) = PP;
qq(1) = 0;
for ind = 1:length(tt)-1
    SS(:,:,ind+1) = SS(:,:,ind)+(AA'*SS(:,:,ind)+SS(:,:,ind)*AA-2*SS(:,:,ind)*BB/R*BB'*SS(:,:,ind)+QQ)*h;
    qq(ind+1) = qq(ind) + trace(SS(:,:,ind)*BB*noise_cov*BB')*h;
end

%% dimension = 2
if d == 2
    
%% Plot result
    
    [gridx, gridy] = meshgrid(gridT{1},gridT{2});
    lqr_values = exp(-(sum(([gridx(:) gridy(:)]*SS(:,:,1)).*[gridx(:) gridy(:)],2)+qq(1))/lambda);
    
    figure
    ht_pdf = surf(gridT{1},gridT{2},reshape(lqr_values,n(2),n(1)),'EdgeColor', 'none');
    colorbar;
    view(0,90)
%     caxis([0 1])
    xlabel('x')
    ylabel('y')
    title(['Time = ',num2str(tt(1))]);
    axis ([bdim(1,:),bdim(2,:)])
    for k = 2:10:ind%length(tt)
         lqr_values = exp(-(sum(([gridx(:) gridy(:)]*SS(:,:,ind)).*[gridx(:) gridy(:)],2)+qq(k))/lambda);
         set(ht_pdf, 'ZData', reshape(lqr_values,n(2),n(1)) );
         title(['Time = ',num2str(tt(k))]);
         drawnow
         pause(1.0/1000);
    end
    
    
    %%
    
    ten_ori = zeros(length(tt),1);
    for k = 1:ind%length(tt)
        ten_ori(k) = EvalT(F_all{k}*F_scale(k),[0,0], gridT);% max(double(F_all{k}));%
    end
    
    figure
    plot(tt,[-lambda*log(ten_ori) qq])
    
    figure
    plot(tt,[ten_ori exp(-qq/lambda)])
    
    figure
    plot(tt,(-lambda*log(ten_ori)-qq))
    
end