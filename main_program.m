% MAIN_PROGRAM - main program of LHJB Toolbox. Solves the linear HJB
% equation.
%
%   INPUT PARAMETERS:
%
%
%   1. DOMAIN
%
%   d - number of state space dimensions
%
%   n - number of grid points in each dimension. Will have the same number
%   of grid points in each dimension if scalar.
%
%   bdim - the dimensions of the hyperrectangle domain; bdim(i,1) and
%   bdim(i,2) is the lower and upper boundary in dimension i.
%
%   bcon - the boundary conditions. bcon{i} is boundary conditions in
%   dimension i:
%   bcon{i} = {'p'}. Periodic
%   bcon{i} = {'v'}. Vanishing boundary
%   bcon{i} = {'d',val_lo,val_up}. Dirichlet with val_lo and val_up 
%   for lower and upper boundary in dimension i.
%   bcon{i} = {'n',val_lo,val_up}. Neumann with val_lo and val_up for
%   lower and upper boundary in dimension i.
%   where all the values are for the desirability function.
%
%   bsca - boundary scaling. bsca(i,1) and bsca(i,2) is how much lower and
%   upper boundary in dimension i should be scaled. The operator is scaled
%   accordingly. Leave empty if scaling method is used.
%
%   region - the dimensions of the goal region given as for bdim.
%
%   regval - the value in the region for the desirability function. Usually
%   1 corresponding to zero cost for the value function.
%
%   regsca - how much the goal region should be scaled. The operator is 
%   scaled accordingly. Leave empty if scaling method is used.
%   
%   sca_ver - specifies the version of the function that creates the
%   boundary scaling if not manually set. sca_ver = 1 is make_bc_sca
%   (second scaling method) and sca_ver = 2 is make_bc_sca_ver (first
%   scaling method).
%
%   ord_of_acc - the order of accuracy for the derivatives. ord_of_acc(i,j)
%   is the order of accuracy for the i:th derivative for interior
%   derivatives if j=1 and edge derivatives if j=2. If ord_of_acc is a
%   scalar, the same order_of_accarucy is set for all derivatives.
%
%
%   2. DYNAMICS
%
%   lambda - lambda in the linear HJB equation.
%
%   x - symblic vector for the state space. Used to write the functions in
%   symbolic notation.
%
%   f - function f from linear HJB equation in symbolic notation.
%   
%   G - function G from linear HJB equation in symbolic notation.
%
%   B - function B from linear HJB equation in symbolic notation.
%
%   noise_cov - the variance vector for the noise corresponding to the
%   diagonal of the covariance matrix (non-diagonal elements are zero).
%   Scalar input is treated as same variance for all dimensions.
%
%   q - the function q from linear HJB equation in symbolic notation.
%
%   R - the control cost matrix from linear HJB equation. R=r*I if input is
%   a scalar r
%
%
%   3. ARTIFICIAL DIFFUSION, COMPRESS OPERATOR AND SOLVE ALS
%   
%   artdiff_option - cell with the following data:
%   artdiff_version - version of the artificial diffusion (AD) term.
%   Possible choices are 1,2 and 3.
%   artdiff_dims - artdiff_dims = 'all' if all dimension shall have AD.
%   artdiff_dims = 'needed' if just dimensions that has a convection term
%   and no corresponding diffusion term shall have AD.
%   artdiff_scale - scale parameter of the AD term.
%
%   tol_err_op - accepted error when compressing operator.
%
%   comp_options - cell with the following data used when compressing
%   operator:
%   tol_it_op - tolerated maximum number of iterations.
%   r_tol_op - minimum decrease in error, adds rank when below.
%   alpha_op - regularization term.
%   if comp_options is empty, default setup will be used.
%
%   tol_err - accepted error when solving system.
%
%   als_options - cell with the following data used when solving system:
%   tol_it - maximum toleranced number of iterations.
%   tol_rank - maximum tolerance rank of solution.
%   e_type - error type, 'average' or 'total' error.
%   r_tol - minimum decrease in error, adds rank when below.
%   alpha - regularization term.
%   newcond_r_tol - tolerated decrease in error during preconditioning.
%   newcond_it_tol - maximum tolerated iterations during preconditioning.
%   if als_options is empty, default setup will be used.
%
%   als_variant - cell with options when running als_sys_var:
%   max_osc - maximum number of oscillations before it start over.
%   max_it_rank - maximum rank for the iterated F.
%   Leave it empty if original als_sys should be run. 
%
%   debugging - 1 for extra documentation from als run (takes extra memory).
%
%
%   4. VISUALIZE RESULTS, RUN SIMULATION AND RUN INDEX
%
%   saveplots - 1 if plots shall be saved.
%
%   savedata - 1 if data shall be saved.
%
%   sim_config - cell with the following data for simulation:
%   T - end time for simulations.
%   x0_list - matrix with starting points; each column corresponds to a
%   starting point.
%   new_noise_cov - variance vector for noise in simulation. The noise
%   from main_run_data will be used if new_noise_cov is empty.
%   new_region - goal region in simulation, replaces the old region. The
%   old region will be used if new_region is empty.
%   no simulation will be run if sim_config is empty.
%
%   5. RUN INDEX
%
%   run - run index, used for multiple runs with main_run.
%
%   OBS:
%   Neumann conditions has not been tested with any example yet.

% Elis Stefansson, Aug 4 2015

% Add the scripts to path - run this once before running any setup
% addpath('./src/');

%% setup

% run('./examples/Inverted_Pendulum.m')
% run('./examples/Smooth_2D_Example.m')
% run('./examples/Quadcopter.m')

%% arrange input
input1 = {d,n,bdim,bcon,bsca,region,regval,regsca,sca_ver,ord_of_acc};
input2 = {x,f,G,B,noise_cov,q,R,lambda};
input3 = {artdiff_option,tol_err_op,comp_options,tol_err,als_options,als_variant,debugging};
input4 = {saveplots,savedata,sim_config};

%% run
[F] = main_run(input1,input2,input3,input4,run);
