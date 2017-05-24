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

nd = length(AtG);

nf = size(AtG{1},1);
rF = ncomponents(F);
rA = round(size(AtA{1},1)/nf);
rG = size(AtG{1},2);

%%% documentation %%%
F_cell = {}; %for no documentation, these variabels remain empty.
B_cell = {};
b_cell = {};
%%% documentation %%%

% Distribute normalization factor evenly
F_U = cell(nd,1);
Flambda = (F.lambda').^(1/nd);
for i = 1:nd
    F_U{i} = bsxfun(@times,F.U{i},Flambda);
end

for k = 1:nd
    idx = [1:(k-1) (k+1):nd];
    nf = size(F_U{k},1);
    B = zeros(rF*nf,rF*nf);
    b = zeros(rF*nf,1);
    
    %Calculate B matrix  
    for i=1:rF      
        for j=1:rF         
            Mt = zeros(nf,nf);
            for ia=1:rA
                for ja=1:rA 
                    Aprod2 = 1;
                    for dd = idx
                        nnf = size(F_U{dd},1);
                        A_ia_ja = AtA{dd}((ia-1)*nnf+(1:nnf),(ja-1)*nnf+(1:nnf));
                        Aprod2 = Aprod2*F_U{dd}(:,j)'*A_ia_ja*F_U{dd}(:,i);
                    end  
                    Mt = Mt + AtA{k}((ia-1)*nf+(1:nf),(ja-1)*nf+(1:nf))'*Aprod2;                   
                end
            end
            
            if i==j
                B((i-1)*nf+(1:nf),(j-1)*nf+(1:nf)) = Mt+alpha*eye(nf);
            else
                B((i-1)*nf+(1:nf),(j-1)*nf+(1:nf)) = Mt;
            end
            
        end
        
        %Calculate b matrix
        temp = ones(1,rG);
        for dd = idx
            temp = temp.*(F_U{dd}(:,i)'*AtG{dd});
        end
        b((i-1)*nf+(1:nf),:) = AtG{k}*temp';
    end
    
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
 
    F_U{k} = reshape(u,nf,rF);
    
    %%% documentation %%%
    if debugging == 1
        F = ktensor(F_U);
        F = arrange(F);
        F_cell{k} = F;
    end
    %%% documentation %%%
    
    
end

F = ktensor(F_U);
F = fixsigns(arrange(F));

warning('on', 'MATLAB:nearlySingularMatrix');
