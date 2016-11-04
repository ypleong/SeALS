%% Example - Smooth 2D Example
% Elis Stefansson, Aug 2015

% The following example is taken from the paper:
% Linearly Solvable Stochastic Control Lyapunov Functions 
% by Yoke Peng Leong et. al.

%% setup

% 1. domain 
d = 2;
x = sym('x',[d,1]); %do not change
n = 101;
bdim = [-1 1 ; -1 1];
bcon = { {'d',exp(-5),exp(-5)} , {'d',exp(-5),exp(-5)} };
bsca = []; %no manual scaling
region = [0 0];
regval = 1;
regsca = [];
sca_ver = 1; %first scaling method
ord_of_acc = 6;

% 2. dynamics
lambda = 1;
f = 2*[x(1)^5-x(1)^3-x(1)+x(1)*x(2)^4; x(2)^5-x(2)^3-x(2)+x(2)*x(1)^4];
G = [x(1) 0; 0 x(2)];
B = G;
noise_cov = 1;
q = x(1)^2+x(2)^2;
R = 1;

% 3. artificial diffusion, compress operator and solve als
artdiff_option = [];
tol_err_op = 1e-6;
comp_options = [];
tol_err = 1e-9;
als_options = [];
als_variant = {10,50};
debugging = 1;

% 4. visualize result and run simulation
saveplots = 1;
savedata = 1;
sim_config = {5,0.005,[-0.5 0.5 0.2 -0.7 -0.2 0.3;0.5 0 0.6 -0.2 -0.4 -0.6],[],[]};

% 5. run index
run = 1;
