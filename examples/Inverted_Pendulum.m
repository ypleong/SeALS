%% Example - Inverted Pendulum
% Elis Stefansson, Aug 2015

% The following example is taken from the paper:
% Linear Hamilton-Jacobi-Bellan Equations in High Dimensions
% by Matanya B. Horowitz et al.

%% setup

% 1. domain 
d = 2;
x = sym('x',[d,1]); %do not change
n = 201;
bdim = [-pi pi ; -11 11];
bcon = { {'p'} , {'d', exp(-50), exp(-50)} };
bsca = []; %no manual scaling
region = 0.3*[-0.6 0.6 ; -1 1];
regval = 1;
regsca = [];
sca_ver = 2; %second scaling method
ord_of_acc = 6;

% 2. dynamics
lambda = 0.2;
m = 2; M = 8; l = .5; g = 9.8; mr = m/(m+M); den = 4/3-mr*cos(x(1))^2;
f1 = (g/l)*sin(x(1))/den;
f2 = -0.5*mr*x(2)^2*sin(2*x(1))/den;
f = [x(2) ; f1+f2];
G = [0 ; -mr/(m*l)*cos(x(1))/den];
B = G;
noise_cov = 10;
q = 0.1*x(1)^2+0.05*x(2)^2;
R = 0.02;

% 3. artificial diffusion, compress operator and solve als
artdiff_option = {1,'needed',1};
tol_err_op = 1e-7;
comp_options = [];
tol_err = 0;
als_options = [];
als_variant = {10,20}; %for speed
debugging = 0;

% 4. visualize result and run simulation
saveplots = 1;
savedata = 1;
sim_config = {60,0.005,repmat([-0.5; 0],1,10),[],[]};

% 5. run index
run = 1;
