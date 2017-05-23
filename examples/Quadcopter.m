%% Example - Quadcopter
% Elis Stefansson, Nov 10 2015

% The following example is taken from the paper:
% Linear Hamilton-Jacobi-Bellan Equations in High Dimensions
% by Matanya B. Horowitz et al.

% The states are:
% [dx, dy, dz, dpsi, dth, dphi, x, y, z, psi, th, phi]
%  1   2   3   4     5    6     7  8  9  10   11  12

%% setup

% 1. domain 
d = 12;
x = sym('x',[d,1]); %do not change
n = 100;
bdim = [-8 8 ; -8 8; -8 8; -10*pi 10*pi; -10*pi 10*pi; -10*pi 10*pi; -1 1; -1 1; -1 1; -pi pi; -pi pi; -pi pi];

% help function for boundary x = 1
bcon_x = (1-x(1)^2/8^2)*(1-x(2)^2/8^2)*(1-x(3)^2/8^2)*(1-x(4)^2/(10*pi)^2)*(1-x(5)^2/(10*pi)^2)*(1-x(6)^2/(10*pi)^2)*(1-x(8)^2)*(1-x(9)^2)*(1-x(10)^2/pi^2)*(1-x(11)^2/pi^2)*(1-x(12)^2/pi^2);

bcon = { {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,0} , {'d',0,bcon_x} , {'d',0,0} , {'d',0,0} , {'p'} , {'p'} , {'p'} };

bsca = []; %no manual scaling
region = []; %no goal region
regval = [];
regsca = [];
sca_ver = 1; %first scaling version
ord_of_acc = 2; %since als error probably is dominant

% 2. dynamics
g = 9.8;
f = [0; 0; -g; 0; 0; 0; x(1); x(2); x(3); x(4); x(5); x(6)];

%help functions for G
hf1 = sin(x(12))*sin(x(10))+cos(x(12))*cos(x(10))*sin(x(11));
hf2 = cos(x(12))*sin(x(11))*sin(x(10))-cos(x(10))*sin(x(12));
hf3 = cos(x(11))*cos(x(12));

G = [hf1 0 0 0; hf2 0 0 0; hf3 0 0 0; ...
    0 1 0 0; 0 0 1 0; 0 0 0 1; ...
    0 0 0 0; 0 0 0 0; 0 0 0 0; ...
    0 0 0 0; 0 0 0 0; 0 0 0 0];

B = G;
noise_cov = 1;
q = 1;
R = 1;
lambda = noise_cov*R;

% 3. artificial diffusion, compress operator and solve als
artdiff_option = {1,'needed',1};
tol_err_op = 1e-5;
comp_options = [];
tol_err = 1e-9;
als_options = [];
als_variant = {10,50};
debugging = 0;

% 4. visualize result, run simulation
saveplots = 1;
savedata = 1;
sim_config = {5,0.005,[0;0;0; 0;0;0; 0;0;0; 0.1;0.1;0.1],[],[]};

% 5. run index
run = 1;