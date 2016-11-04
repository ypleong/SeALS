%% Example - Scalar Unstable System
% Elis Stefansson, Aug 2015

% The following example is taken from the paper:
% Linearly Solvable Stochastic Control Lyapunov Functions 
% by Yoke Peng Leong et. al.

%% setup

% 1. domain 
d = 1;
n = 101;
bdim = [-1 1];
bcon = { {'d',20*exp(-10),20*exp(-10)} };
bsca = []; %no manual scaling
region = [0 0];
regval = 1;
regsca = [];
sca_ver = 2; %second scaling method
ord_of_acc = 6;

% 2. dynamics
lambda = 1;
x = sym('x',[d,1]); %do not change
f = x(1)^3+5*x(1)^2+x(1);
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
debugging = 1;

% 4. visualize result and run simulation
saveplots = 1;
savedata = 1;
sim_config = {5,0.005,[-0.5 0.5],[],[]};

% 5. run index
run = 1;

