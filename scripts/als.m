function [F, err, iter, Fcond, e_list, t_step, illcond, noreduce] = als(G,F,e,options)
% ALS code to perform the ALS algorithm as described in Beylkin
% and Mohlenkamp 2005.
% Input:
%   G - tensor to be approximated.
%   F - inital guess, if empty a random guess will be used.
%   e - desired accuracy.
%   options - cell with the following data:
%      maxit - maximum number of iterations.
%      r_tol - minimum decrease for error, adds rank when below.
%      alpha - regularization term.
%      if options is empty, ALS will use the default setup maxit = 2000, r_tol
%      = 1e-3 and alpha = 1e-14.
% Output:
%   F - the low rank approimxation of G.
%   err - the error during the run.
%   iter - the number of iterations.
%   Fcond - the condition number of F during the run.
%   e_list - index when tensor term is added.
%   t_step - computing time for als onestep during the run.
%   illcond - 1 if a als matrice were illcond, otherwise 0.
%   noreduce - 1 if als couldn't reduce rank, otherwise 0.
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 2015
% new version for Toolbox

if isempty(options) == 1
    maxit = 2000;
    r_tol = 1e-3; %when to add to the rank
    alpha = 1e-14; %regularization
else
    [maxit,r_tol,alpha] = deal(options{:});
end

noreduce = 0;
illcond = 0;

G = arrange(G);
nd = ndims(G);
sizeG = size(G);
rG = ncomponents(G);
Fcond = 1;

if isempty(F)
    U = cell(1,nd);
    for n = 1:nd
        U{n} = randn(sizeG(n),1);
        U{n} = U{n}/norm(U{n});
    end
    F = ktensor(U);
end
rF = ncomponents(F);

Flambda = F.lambda;
Glambda = G.lambda;

FU = cell(nd,1);
GU = cell(nd,1);
for k = 1:nd
    FU{k} = F.U{k};
    GU{k} = G.U{k};
end

norG = norm(G);
old_err = norm(F-G)/norG;
%old_err = norm(F-G);

reverseStr = '';

e_count = 2;
e_list(1) = 1;

for iter = 1:maxit
    
    %display progress
    if 1 == 0
        msg = sprintf('Iteration: %d (%d), error=%2.3f (%2.3f)', iter, maxit, old_err, e); %don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    step_time = tic;
    [Flambda, FU, status] = als_onestep(Glambda,GU,Flambda,FU,alpha);
    t_step(iter) = toc(step_time);
    
    if status == 0
        fprintf('ALS matrix inversion ill-conditioned, returning after iteration %d\n', iter);
        illcond = 1;
        break;
    end
    
    F = arrange(F);
    Fcond(iter) = norm(F.lambda)/norm(F);
    
    ETlambda = [Flambda ; -1*Glambda];
    ETU = cell(nd,1);
    for k = 1:nd
        ETU{k} = [FU{k} GU{k}];
    end
    ET = ktensor(ETlambda,ETU);
    
    err(iter) = norm(ET)/norG;
    %err(iter) = norm(ET);
    
    if err(iter) <= e
        break;
    end
    
    if abs(err(iter) - old_err)/old_err < r_tol
        clear nF U
        rF = rF + 1;
        e_list(e_count) = iter;
        e_count = e_count + 1;
        
        if rF == rG
            F = G;
            err(iter) = 0;
            disp('cannot reduce rank')
            noreduce = 1;
            break;
        end
        
        for n = 1:nd
            U = randn(sizeG(n),1);
            U = U/norm(U);
            FU{n} = [FU{n} U];
        end
        Flambda = [Flambda; 1];
        
    end
    old_err = err(iter);
end

if (iter == maxit && err(iter) > e)
    disp('Max iterations reached without finding a good solution')
    noreduce = 1;
    F = G;
end

if ~noreduce
    F = ktensor(Flambda,FU);
end