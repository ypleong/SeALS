function [M_cell,Aprodsum_cell] = als_mat_est(A,F)
% ALS_MAT_EST estimates the matrices M{1,1} and the produkt terms in
% M{1,1} for each dimension, which is used to scale the boundary conditions
% in MAKE_BC_SCA. ALS_MAT_EST originates from ALS_ONESTEP_SYS. Thus, see
% ALS_ONESTEP_SYS first. See also the paper.
% Inputs:
%   A - the operator.
%   F - estimated solution F to AF=G.
% Outputs:
%   M_cell - cell with all the matrices M{1,1} in each dimension.
%   Aprodsum_cell - cell with the products terms in M{1,1} in each
%   dimension.
%
% See also MAKE_BC_SCA.

% Elis Stefansson, Aug 2015

A = arrange(A);
F = arrange(F);

nd = ndims(F);

sizeA = size(A);
rA = ncomponents(A);
rF = ncomponents(F);

% incooporate normalization factor, evenly distributed;
for i = 1:nd
    A.U{i} = A.U{i}*diag((A.lambda).^(1/nd));
    F.U{i} = F.U{i}*diag((F.lambda).^(1/nd));
end
A.lambda = ones(rA,1);
F.lambda = ones(rF,1);

% get matrices
for k = 1:nd
    nf = sqrt(sizeA(k));
    Acell{k} = reshape(A.U{k},nf,nf,rA);
end

%%% for scaling %%%
M_cell = cell(1,nd);
Aprodsum_cell = cell(1,nd);
%%% for scaling %%%

for k = 1:nd
    nf = sqrt(sizeA(k));
    
    idx = [1:(k-1) (k+1):nd];
    Fcell = cell(nd,1);
    for n = 1:nd
        Fcell{n} = F.U{n};
    end
    
    M = cell(rF,rF);
    
    %Calculate B matrix
    for i=1:rF
        for j=1:rF
            
            M{i,j} = zeros(nf);
            
            %%% for scaling %%%
            Aprodsum = 0;
            %%% for scaling %%%
            
            %%% for scaling %%%%
            % ignored for scaling case
            %if i==j
            %    M{i,j} = M{i,j}+alpha*eye(nf);
            %end
            %%% for scaling %%%
            
            for ia=1:rA
                A_k_ia = Acell{k}(:,:,ia);
                for ja=1:rA
                    A_k_ja = Acell{k}(:,:,ja);
                    
                    Aprod = 1;
                    for d = idx
                        
                        A_d_ia = Acell{d}(:,:,ia);
                        A_d_ja = Acell{d}(:,:,ja);
                        
                        u_d_i = Fcell{d}(:,i);
                        u_d_j = Fcell{d}(:,j);
                        
                        Aprod = Aprod * (A_d_ia*u_d_j)' * (A_d_ja*u_d_i);
                    end
                    
                    Aprodsum = Aprodsum+Aprod;
                    M{i,j} = M{i,j} + A_k_ja'*A_k_ia*Aprod;
                    
                end
            end
            
        end
    end
    
    %%% for scaling %%% 
    M_cell{k} = M{1,1};
    Aprodsum_cell{k} = Aprodsum;
    %%% for scaling %%%
    
end