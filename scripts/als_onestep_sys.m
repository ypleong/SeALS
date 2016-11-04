function [F, status, F_cell, B_cell, b_cell] = als_onestep_sys(G,F,A,alpha,debugging)
% ALS_ONSTEP_SYS performs one step of the ALS algorithm for solving AF=G.
% Based on the version in Baylkin and Mholenkamp 2005.
% Inputs:
%   G - the RHS as a ktensor
%   F - the current solution as a ktensor.
%   A - the operator as a ktensor
%   alpha - regularization constant.
%   debugging - 1 if all F, B and b are saved in cell arrays.
% Outputs:
%   F - updated F.
%   status - 0 if ALS matrices became ill-conditioned. Otherwise 1.
%   F_cell - all F so far during the run.
%   B_cell - all B (ALS matrices) so far during the run.
%   b_cell - all b (RHS vectors) so far during the run.
%
% See also ALS_SYS, ALS_SYS_VAR.

% Elis Stefansson, Aug 2015

status = 1;

G = arrange(G);
F = arrange(F);
A = arrange(A);

nd = ndims(G);
sizeG = size(G);
rG = ncomponents(G);
rF = ncomponents(F);
rA = ncomponents(A);

%%% documentation %%%
F_cell = {}; %for no documentation, these variabels remain empty.
B_cell = {};
b_cell = {};
%%% documentation %%%

%This let's us not worry about the normalization factor. It might not be a
%good idea in term of condition to put all the scaling on the first
%element.
% A.U{1} = A.U{1}*diag(A.lambda);
% G.U{1} = G.U{1}*diag(G.lambda);
% F.U{1} = F.U{1}*diag(F.lambda);

%Let's distribute it evenly, see if it makes a difference;
for i = 1:nd
    A.U{i} = A.U{i}*diag((A.lambda).^(1/nd));
    G.U{i} = G.U{i}*diag((G.lambda).^(1/nd));
    F.U{i} = F.U{i}*diag((F.lambda).^(1/nd));
end

A.lambda = ones(rA,1);
G.lambda = ones(rG,1);
F.lambda = ones(rF,1);

% tensor reference is slow, so do this at the expense of mem
for k = 1:nd
    nf = sizeG(k);
    Acell{k} = reshape(A.U{k},nf,nf,rA);
    Gcell{k} = G.U{k};
end

for k = 1:nd
    nf = sizeG(k);
    
    idx = [1:(k-1) (k+1):nd];
    Fcell = cell(nd,1);
    for n = 1:nd
        Fcell{n} = F.U{n};
    end
    
    M = cell(rF,rF);
    N = cell(rF,1);
    
    %Calculate B matrix
    for i=1:rF
        for j=1:rF
            
            M{i,j} = zeros(nf);
            
            if i==j
                M{i,j} = M{i,j}+alpha*eye(nf);
            end
            
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
                    
                    M{i,j} = M{i,j} + A_k_ja'*A_k_ia*Aprod;
                    
                end
            end
            
        end
    end
    
    %Calculate b matrix
    for i=1:rF
        N{i} = zeros(nf,1);
        
        for ia=1:rA
            A_k_ia = Acell{k}(:,:,ia);
            for ig=1:rG
                g_k_ig = Gcell{k}(:,ig);
                Agprod = 1;
                for d=idx%setdiff(1:nd,k)
                    A_d_ia = Acell{d}(:,:,ia);
                    
                    g_d_ig = Gcell{d}(:,ig);
                    
                    u_d_i = Fcell{d}(:,i);
                    
                    Agprod = Agprod * (A_d_ia*u_d_i)'*(g_d_ig);
                end
                
                N{i} = N{i} + A_k_ia'*g_k_ig*Agprod;
            end
        end
    end
    
    B = cell2mat(M);
    b = cell2mat(N);
    
    %%% Debugging %%%
    if 1==0 %"movie-plot" of ALS matrices
        plotmatrix(B);
        %error('done')
    end
    %%% Debugging %%%
    
    %%% documentation %%%
    if debugging == 1
        B_cell{k} = B;
        b_cell{k} = b;
    end
    %%% documentation %%%
    
    if cond(B) > 1e13
        status = 0;
        %%% Debugging %%%
        if 1==0
            plotmatrix(B);
        end
        %%% Debugging %%%
    end
    
    %disp(B);
    u = B\b;
    
    u = reshape(u,nf,rF);
    newU = cell(d,1);
    
    for i=1:d
        newU{i} = F.U{i};
    end
    
    newU{k} = u;
    
    F = ktensor(newU);
    
    %%% documentation %%%
    if debugging == 1
        F_cell{k} = F;
    end
    %%% documentation %%%
    
    F = arrange(F);
    
    % again, so we don't have to worry about normalization
    for i = 1:nd
        F.U{i} = F.U{i}*diag((F.lambda).^(1/nd));
    end
    F.lambda = ones(rF,1);
    
end

F = arrange(F);