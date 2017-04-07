function [F, status, F_cell, B_cell, b_cell] = als_onestep_sys(AtA,AtG,F,alpha,debugging)
% ALS_ONSTEP_SYS performs one step of the ALS algorithm for solving AF=G.
% Based on the version in Baylkin and Mholenkamp 2005.
% Inputs:
%   AtG - A'*G (nf,rA*rG,nd)
%   F - the current solution as a ktensor.
%   AtA - A'*A (nf*rA,nf*rA,nd)
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

warning('error', 'MATLAB:nearlySingularMatrix');

status = 1;

F = arrange(F);

sizeAtG = size(AtG);
if length(sizeAtG) == 3
    nd = sizeAtG(3);
elseif length(sizeAtG) == 2
    nd = 1;
end
nf = sizeAtG(1);

rF = ncomponents(F);
rA = round(size(AtA,1)/nf);

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
F_U = zeros(nf,rF,nd);
for i = 1:nd
    F_U(:,:,i) = F.U{i}*diag((F.lambda).^(1/nd));
end

for k = 1:nd
    idx = [1:(k-1) (k+1):nd];
    
    M = cell(rF,rF);
    N = cell(rF,1);
    
    %Calculate B matrix    
    for i=1:rF
        for j=1:rF
            
            M{i,j} = zeros(nf);
            
            if i==j
                M{i,j} = M{i,j}+alpha*eye(nf);
            end
            
%             AtA2 = permute(reshape(AtA,nf,rA,nf,rA,nd),[ 2 4 3 1 5]);
%             tic
%             F_i2 = permute(repmat(F_U(:,i,:),1,1,1,rA,rA),[4 5 2 1 3]);
%             F_j2 = permute(repmat(F_U(:,j,:),1,nf,1,rA,rA),[4 5 1 2 3]);
%             temp2 = squeeze(dot(dot(F_j2, AtA2,3),F_i2,4));
%             MM = repmat(prod(temp2(:,:,idx),3),1,1,nf,nf);
%             M2 = M{i,j}+squeeze(sum(sum(MM.* permute(AtA2(:,:,:,:,k),[1 2 4 3 5]),1),2));
%             toc
            
%             tic
            F_i = squeeze(F_U(:,i,:));
            F_j = repmat(F_U(:,j,:),1,nf,1);      
            for ia=1:rA
                for ja=1:rA                    
                    A_ia_ja = AtA((ia-1)*nf+(1:nf),(ja-1)*nf+(1:nf),:);
                    temp = dot(squeeze(dot(F_j, A_ia_ja)), F_i);               
                    Aprod2 = prod(temp(idx));
                    
                    M{i,j} = M{i,j} + A_ia_ja(:,:,k)'*Aprod2;
                    
                end
            end
%             toc
            
        end
    end
    
    %Calculate b matrix
    for i=1:rF
        temp = dot(repmat(F_U(:,i,:),1,sizeAtG(2),1), AtG);
        N{i} = sum(AtG(:,:,k).*repmat(prod(temp(:,:,idx),3),nf,1),2);
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
    
    try
        u = B\b;
    catch
        status = 0;
        
        %%% Debugging %%%
        if 1==0
            plotmatrix(B);
        end
        %%% Debugging %%%
        warning('off', 'MATLAB:nearlySingularMatrix');
        u = B\b;
        warning('error', 'MATLAB:nearlySingularMatrix');
    end
    
    u = reshape(u,nf,rF);
    
    F_U(:,:,k) = u;
    
    %%% documentation %%%
    if debugging == 1
        U = cell(nd,1);
        for ii = 1:nd
            U{ii} = F_U(:,:,ii);
        end
        F = ktensor(U);
        F = arrange(F);
        F_cell{k} = F;
    end
    %%% documentation %%%
    
    
end

U = cell(nd,1);
for ii = 1:nd
    U{ii} = F_U(:,:,ii);
end
F = ktensor(U);
F = arrange(F);

warning('on', 'MATLAB:nearlySingularMatrix');
