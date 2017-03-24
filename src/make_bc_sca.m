function [bsca, regsca] = make_bc_sca(op,bcon,region,regval,als_options,fd1,grid,x,n)
% MAKE_BC_SCA creates scalings for the boundary condition and goal region
% using estimations of unscaled matrices in ALS_SYS. See also the paper.
% Input:
%   op - the operator as ktensor.
%   bcon - boundary conditions specifications, see MAIN_PROGRAM for
%   details.
%   region - the dimension of the goal region.
%   als_options - cell withoptions for the als run, see MAIN_PROGRAM for 
%   details.
%   fd1 - cell of first order differentiation matrices. fd1{i} is
%   differentiation in dimension i.
%   grid - the discretization grid.
%   x - the symbolic state space vector.
%   n - number of grid points in each dimension.
% Outputs:
%   bsca - boundary scaling. bsca(i,1) and bsca(i,2) is how much lower and
%   upper boundary in dimension i should be scaled.
%   regsca - how much the goal region should be scaled.
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 2015

d = length(n);

%% Step 1 - create unscaled op and bc

% unscaled
bsca = ones(d,2);
regsca = 1;

% bc
[bc] = makebc(bcon,bsca,grid,x,n);

% operator
[op] = makebcop(op,bcon,bsca,n,fd1);
if isempty(region) == 0
    [op,bc] = incorpregion(op,bc,region,grid,regval,regsca);
end

%% Step 2 - estimate rank 1 solution

if isempty(als_options) == 1
    als_options = {2000,1,'average',1e-3,1e-12,0.01,15};
else
    als_options{2} = 1;
end

fprintf('Estimating rank 1 solution for boundary scaling. \n')
F_est = als_sys(op,bc,[],0,als_options,0);
fprintf('Estimation complete. \n')

%% Step 3 - get scaling from estimated als matrices

% scaling values
bsca = ones(d,2);
regsca = 0;

% obtain relevant als matrices
[M_cell,Aprodsum_cell] = als_mat_est(op,F_est);

% get scale factors for each dim
for i = 1:d
    
    bconi = bcon{i};
    
    % get block matrix M{1,1} and its product sum
    M_11 = M_cell{i};
    Aprodsum = Aprodsum_cell{i};
    
    % get diagonals
    du = abs(diag(M_11,1));
    dd = abs(diag(M_11,0));
    dl = abs(diag(M_11,-1));
    
    % for lower boundary
    sca_lo = sum( du(2:3)+dd(2:3)+dl(2:3) )/6; %average magnitud
    sca_lo = sca_lo/Aprodsum; %divide by how much matrices in M_11 has been scaled
    sca_lo = sqrt(sca_lo); %since the summation of M_11 has two matrices multiplied
    
    % for upper boundary
    sca_up = sum( du(end-2:end-1)+dd(end-2:end-1)+dl(end-2:end-1) );
    sca_up = sca_up/Aprodsum;
    sca_up = sqrt(sca_up);
    
    if strcmp(bconi{1},'p') == 1 %then must have the same scalings
        bsca(i,1) = (sca_lo+sca_up)/2;
        bsca(i,2) = (sca_lo+sca_up)/2;
    else
        bsca(i,1) = sca_lo;
        bsca(i,2) = sca_up;
    end
    
    % for region
    if isempty(region) == 0
        gridi = grid{i};
        
        reg_idx = find(gridi>=region(i,1) & gridi<=region(i,2));
        reg_min = min(reg_idx);
        reg_max = max(reg_idx);
        
        s1 = sum( du(reg_min-2:reg_min-1)+dd(reg_min-2:reg_min-1)+dl(reg_min-2:reg_min-1) );
        s2 = sum( du(reg_max+1:reg_max+2)+dd(reg_max+1:reg_max+2)+dl(reg_max+1:reg_max+2) );
        
        sca_reg = (s1+s2)/12; %take mean
        sca_reg = sca_reg/Aprodsum;
        sca_reg = sqrt(sca_reg);
        
        regsca = regsca+sca_reg;
    end
    
end

if isempty(region) == 0
    regsca = regsca/d;
end








