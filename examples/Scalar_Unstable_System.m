%% Example - Scalar Unstable System
% Elis Stefansson, Aug 2015

% The following example is taken from the paper:
% Linearly Solvable Stochastic Control Lyapunov Functions 
% by Yoke Peng Leong et. al.

%% setup

% 1. domain 
d = 1;
n = 151;
bdim = [-5 5];
bcon = { {'d',1,0} };
bsca = []; %no manual scaling
region = [];
regval = 1;
regsca = 1;
sca_ver = 1; %second scaling method
ord_of_acc = 4;

% 2. dynamics
lambda = 1;
x = sym('x',[d,1]); %do not change
f = x(1);%x(1)^3+5*x(1)^2+x(1);
G = 1;
B = 1;
noise_cov = 1;
q = x(1)^2;
R = 1;

% 3. artificial diffusion, compress operator and solve als
artdiff_option = [];
tol_err_op = 1e-7;
comp_options = [];
tol_err = 10^-6;
als_options = [];
als_variant = []; %org. ALS (SeALS gives the same result bcs AF=G is just a matrix eq. in 1D)
debugging = 0;

% 4. visualize result and run simulation
plotdata = 1;
saveplots = 0;
savedata = 1;
sim_config = {5,0.005,[-0.5 0.5],[],[]};

% 5. run index
run = 1;

% Setup run
input1 = {d,n,bdim,bcon,bsca,region,regval,regsca,sca_ver,ord_of_acc};
input2 = {x,f,G,B,noise_cov,q,R,lambda};
input3 = {artdiff_option,tol_err_op,comp_options,tol_err,als_options,als_variant,debugging};
input4 = {saveplots,savedata,plotdata,sim_config};

save(['./test_run/scalar_',num2str(run),'_setup'])

% Solve for F
[F, grid] = main_run_spectral(input1,input2,input3,input4,run,'./test_run/scalar_');
