function [F, err, iter, Fcond, e_list, time_step, illcond, maxit, maxrank, F_cell, B_cell, b_cell] = als_sys(A,G,F,e,als_options,debugging,verbose)
% ALS_SYS code to perform the ALS algorithm as described in Beylkin
% and Mohlenkamp 2005. Tries to solve AF = G.
% Inputs:
%   F - inital guess. If it is empty a random guess will be used. If
%   it's a scalar ALS_SYS will start with a random tensor of that
%   rank.
%   e - desired error
%   debugging - if 1 all F, B and b are saved in cell arrays.
%   tol_it - maximum number of tolerated iterations. Default if 1000.
%   tol_rank - maximum tolerated rank for F. Default is 20.
%   e_type - type of error ('average' or 'total'). Default is (point)
%   average. 
%   r_tol - minimum decrease in error between iterations. Default is 1e-3.
%   alpha - regularization factor. Default is 1e-12.
%   newcond_r_tol - tolerated decrease for preconditon. Default is 1e-2.
%   newcond_it_tol - maximum number of tolerated iterations for
%   preconditon. Default is 15.
% Outputs:
%   F - achieved solution.
%   err - achieved error.
%   iter - number of iterations.
%   Fcond - condition number for F during the run.
%   e_list - indices when new tensor term was added
%   t_step - computing time for als onestep during the run
%   illcond - 1 if an als matrix were illcond.
%   maxit - 1 if maximum number of tolerated iterations was exceeded.
%   F_cell - all F during a run
%   B_cell - all B (ALS matrices) during a run
%   b_cell - all b (RHS vectors) during a run
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 2015

%% assign
if isempty(als_options)
    tol_it = 2000;
    tol_rank = 30;
    e_type = 'average';
    r_tol = 1e-3;
    alpha = 1e-12;
    newcond_r_tol = 0.01;
    newcond_it_tol = 15;
else
    [tol_it,tol_rank,e_type,r_tol,alpha,newcond_r_tol,newcond_it_tol] = ...
        deal(als_options{:});
end

if strcmp(e_type,'average')
    norA = sqrt(prod(size(G)));
elseif strcmp(e_type,'total')
    norA = 1;
else
    error('wrong type of error specified');
end

if nargin < 7
    verbose = 1;
end

%% documentation
if debugging == 1
    maxit = 0; %1 if max iter exceeded
    maxrank = 0; %1 if max rank exceeded
    illcond = 0; %1 if ill-conditioned matrix
    B_cell{1,2} = {}; %just to make sure it has correct dims
    b_cell{1,2} = {};
else
    maxit = 0;
    maxrank = 0;
    illcond = 0;
    F_cell = {};
    B_cell = {};
    b_cell = {};
end

%% main script
nd = ndims(G);
sizeG = size(G);

U = cell(nd,1);

if isempty(F)
    for n = 1:nd
        U{n} = matrandnorm(sizeG(n),1);
    end
    F = ktensor(U);
    
elseif isfloat(F)
    for n = 1:nd
        U{n} = matrandnorm(sizeG(n),F);
    end
    F = ktensor(U);
else
    F = arrange(F);
    F = fixsigns(F);
end

rF = ncomponents(F);

%%% documentation %%%
if debugging == 1
    F_cell{1,4} = {F};
end
%%% documentation %%%

old_err = norm(SRMultV(A,F)-G)/norA;

reverseStr = '';
msg = '';

useStop = 1;
e_count = 2;
e_list(1) = 1;

if useStop
    FS = stoploop({'Exit execution at current iteration?', 'Stop'}) ;
end

[AtA, AtG] = prepareAG_4_als_sys(A, G);

for iter = 1:tol_it
    
    %%%%%% Debugging %%%%%%
    if length(size(F)) == 0 %set to 2 for 2D-solution plots while solving
        sol = double(F)';
        %[xx,yy] = meshgrid(x1,x2);
        %surf(xx,yy,sol,'EdgeColor','none');
        figure(6)
        %surf(sol,'EdgeColor','none');
        imagesc(sol);
        xlabel('x');ylabel('y');zlabel('u');
        if iter ~= 1
            title(['Solution, time step:', num2str(time_step(iter-1))]);
        end
        pause(0.0001)
    end
    %%%%%% Debugging %%%%%%
    
    % Display progress
    if verbose
        msg = sprintf('Iteration: %d (%d), error=%2.9f (%2.9f), rank(F)=%d\n', iter, tol_it, old_err, e, ncomponents(F));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    step_time = tic;
    [Fn, status, F_cell_onestep, B_cell_onestep, b_cell_onestep] = als_onestep_sys(AtA,AtG,F,alpha,debugging);
    time_step(iter) = toc(step_time);
    
    %%% documentation %%%
    if debugging == 1
        F_cell{iter,1} = F_cell_onestep;
        B_cell{iter,1} = B_cell_onestep;
        b_cell{iter,1} = b_cell_onestep;
    end
    %%% documentation %%%
    
    if status == 0
        fprintf('ALS matrix inversion ill-conditioned, returning after iteration %d\n', iter);
        illcond = 1;
        Fcond(iter) = norm(F.lambda)/norm(F);
        err(iter) = norm(SRMultV(A,F)-G)/norA;
        break;
    else
        F = Fn;
    end
    
    Fcond(iter) = norm(F.lambda)/norm(F);
    err(iter) = norm(SRMultV(A,F)-G)/norA;
    %err(iter) = norm(SRMultV(A,F)-G);
    
    if err(iter) <= e
        break;
    end
    
    if abs(err(iter) - old_err)/old_err < r_tol
        
        rF = rF + 1;
        
        if( rF > tol_rank )
            disp('Maximum rank reached. Quitting...');
            maxrank = 1;
            break;
        end
        
        e_list(e_count) = iter;
        e_count = e_count + 1;
        
        % debugging:
        %fprintf('increased rank on iteration %i\n',iter);
        addStr = repmat(sprintf(' '), 1, length(msg));
        fprintf(addStr);
        
        for n = 1:nd
            U{n} = matrandnorm(sizeG(n),1);
        end
        nF = ktensor(U);
        
        % Precondition the new rank 1 tensor
        count = 1;
        err_newF = 1;
        err_oldF = 2;
        
        % debugging:
        %fprintf('Conditioning new term. Error: %f\n', norm(G-SRMultV(A,F)));
        
        [AtA2, AtG2] = prepareAG_4_als_sys(A, G-SRMultV(A,F));
        
        while count < newcond_it_tol && (abs(err_newF - err_oldF)/err_oldF) > newcond_r_tol
            
            [nF, ~, F_cell_onestep, B_cell_onestep, b_cell_onestep] = als_onestep_sys(AtA2,AtG2,nF,alpha,debugging);

            err_oldF = err_newF;
            err_newF = nF.lambda;
            
            %%% documentation %%%
            if debugging == 1
                if count == 1
                    F_cell_precond = F_cell_onestep;
                    B_cell_precond = B_cell_onestep;
                    b_cell_precond = b_cell_onestep;
                else
                    F_cell_precond = [F_cell_precond, F_cell_onestep];
                    B_cell_precond = [B_cell_precond, B_cell_onestep];
                    b_cell_precond = [b_cell_precond, b_cell_onestep];
                end
            end
            %%% documentation %%%
            
            count = count + 1;
            
        end
        
        %%% documenation %%%
        if debugging == 1
            F_cell{iter,2} = F_cell_precond;
            B_cell{iter,2} = B_cell_precond;
            b_cell{iter,2} = b_cell_precond;
        end
        %%% documentation %%%
        
        F = arrange(F + nF);
        F = fixsigns(F);
        
        %%% documentation %%%
        if debugging == 1
            F_cell{iter,3} = {F};
        end
        %%% documentation %%%
        
        % debugging %
        %fprintf('Conditioning complete. Error: %f, count:%d\n', norm(G-SRMultV(A,F)), count);
        
    end
    old_err = err(iter);
    
    %Check quit option
    if useStop
        if FS.Stop()
            break;
        end
    end
end

if iter == tol_it
    maxit = 1;
end

if useStop
    FS.Clear() ;  % Clear up the box
    clear FS ;    % this structure has no use anymore
end

if ~illcond && verbose
    msg = sprintf('Iteration: %d (%d), error=%2.9f (%2.9f), rank(F)=%d\n', iter, tol_it, err(iter), e, ncomponents(F));
    fprintf([reverseStr, msg]);
end

end