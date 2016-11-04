%% Example - VTOL Aircraft
% Elis Stefansson, Aug 2015

% The following example is taken from the paper:
% Linear Hamilton-Jacobi-Bellan Equations in High Dimensions
% by Matanya B. Horowitz et al.

% The states are:
% [x, dx, y, dy, th, dth]
%  1  2   3  4   5   6

%% setup

% 1. domain 
d = 6;
x = sym('x',[d,1]); %do not change
n = 100;
bdim = [-4 4 ; -8 8 ; 0 2; -1 1; -pi pi; -5 5];

% help function for boundary y = 0
bcon_y = (1-x(1)^2/4^2)*(1-x(2)^2/8^2)*(1-x(4)^2)*(1-x(5)^2/pi^2)*(1-x(6)^2/5^2);

bcon = { {'d',0,0} , {'d',0,0}, {'d',bcon_y,0}, {'d',0,0}, {'p'}, {'d',0,0} };

bsca = [];
region = []; %no goal region
regval = [];
regsca = [];
sca_ver = 1;
ord_of_acc = 2; %since als error is probably dominant

% 2. dynamics
lambda = 6;
g = 9.8; eps = 0.01;
f = [x(2); 0; x(4); -g; x(6); 0];
G = [0 0; -sin(x(5)) eps*cos(x(5)); 0 0; cos(x(5)) eps*sin(x(5)); 0 0; 0 1];
B = G;
noise_cov = 3;
q = 1;
R = 2;

% 3. artificial diffusion, compress operator and solve als
artdiff_option = {1,'needed',0.03};
tol_err_op = 1e-5;
comp_options = [];
tol_err = 1e-9;
als_options = [];
als_variant = {10,50};
debugging = 0;

% 4. visualize result, run simulation and run index.
saveplots = 1;
savedata = 1;
sim_config = {5,0.005,[2; 0; 1; 0; 0.1; 0],[],[]};

%5. run index
run = 1;
