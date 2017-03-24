function [bsca, regsca] = make_bc_sca_var(A,grid,region,bcon)
% MAKE_BC_SCA_VAR creates scalings for the boundary condition and goal
% region using esimations of the matrices from the ktensor of the operator
% op. op is here denoted by A since the code is a modification of
% the code ALS_ONSTEP_SYS.
% Inputs:
%   A - the operator (op).
%   grid - the discratization grid.
%   region - the dimensions of the goal region.
%   bcon - the boundary conditions. bcon{i} is boundary conditions in
%   dimension i:
%      bcon{i} = {'p'}. Periodic
%      bcon{i} = {'d',val_lo,val_up}. Dirichlet with val_lo and val_up 
%      for lower and upper boundary in dimension i.
%      bcon{i} = {'n',val_lo,val_up}. Neumann with val_lo and val_up for
%      lower and upper boundary in dimension i.
%      where all the values are for the desirability function.
%   bsca - boundary scaling. bsca(i,1) and bsca(i,2) is how much lower and
%   upper boundary in dimension i should be scaled.
% Outputs:
%   bsca - boundary scaling. bsca(i,1) and bsca(i,2) is how much lower and
%   upper boundary in dimension i should be scaled.
%   regsca - how much the goal region should be scaled.
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 2015

A = arrange(A);
d = ndims(A);
rA = ncomponents(A);
sizeA = size(A);

% distribute the normalization constants evenly, so one does not need to
% worry about normailization constants.
for i = 1:d
    A.U{i} = A.U{i}*diag((A.lambda).^(1/d));
end
A.lambda = ones(rA,1);

% convert to matrices
for k = 1:d
    nf = sqrt(sizeA(k));
    Acell{k} = reshape(A.U{k},nf,nf,rA);
end

regsca = 0;

for k = 1:d
    
    bconi = bcon{k};
    
    nf = sqrt(sizeA(k));
    matrix_sum = zeros(nf,nf);
    weight_sum = 0;
    
    % Step 1. get weighted matrix sum in dimension i 
    for ia = 1:rA
        
        A_k_ia = Acell{k}(:,:,ia);
        weight = norm(A_k_ia); %scale term with how big it is
        
        matrix_sum = matrix_sum+A_k_ia*weight;
        weight_sum = weight_sum+weight;
        
    end
    matrix_sca = matrix_sum/weight_sum; %take weigthed average
    
    % Step 2. get scalings from matrix_sca
    
    du = abs(diag(matrix_sca,1));
    dd = abs(diag(matrix_sca,0));
    dl = abs(diag(matrix_sca,-1));
    
    % for outer bc
    sca_lo = sum( du(2:3)+dd(2:3)+dl(2:3) )/6; %average
    sca_up = sum( du(end-2:end-1)+dd(end-2:end-1)+dl(end-2:end-1) )/6;
    
    if strcmp(bconi{1},'p') == 1 %then must have the same scalings
        % take mean
        bsca(k,1) = (sca_lo+sca_up)/2;
        bsca(k,2) = (sca_lo+sca_up)/2;
    else
        bsca(k,1) = sca_lo;
        bsca(k,2) = sca_up;
    end
    
    % for region
    if isempty(region) == 0
        gridi = grid{k};
        
        reg_idx = find(gridi>=region(k,1) & gridi<=region(k,2));
        reg_min = min(reg_idx);
        reg_max = max(reg_idx);
        
        s1 = sum( du(reg_min-2:reg_min-1)+dd(reg_min-2:reg_min-1)+dl(reg_min-2:reg_min-1) );
        s2 = sum( du(reg_max+1:reg_max+2)+dd(reg_max+1:reg_max+2)+dl(reg_max+1:reg_max+2) );
        
        sca_reg = (s1+s2)/12; %take mean
        regsca = regsca+sca_reg;
    end
    
end

if isempty(region) == 0
    regsca = regsca/d;
end










