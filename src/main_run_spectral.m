function [F, grid] = main_run_spectral(input1,input2,input3,input4,run,save_folder)
% MAIN_RUN obtains the solution to the setup specified in MAIN_PROGRAM,
% visualizes the result and runs simulations with the obtained controller.
% Inputs:
%   The inputs: input1, input2, input3 and input4 and run corresponds to
%   setup 1,2,3,4 and 5 respectively in MAIN_PROGRAM. See MAIN_PROGRAM for
%   more info.
% Outputs:
%   F - obtained desirability function as ktensor.
%
% See also MAIN_PROGRAM.

% Elis Stefansson, Aug 5 2015

start_whole = tic;

%% Step 1: extract input data

[d,n,bdim,bcon,bsca,region,regval,regsca,sca_ver,ord_of_acc] = deal(input1{:});

[x,f,G,B,noise_cov,q,R,lambda] = deal(input2{:});

[artdiff_option,tol_err_op,comp_options,tol_err,als_options,als_variant,debugging] = deal(input3{:});

[saveplots,savedata,plotdata,sim_config] = deal(input4{:});

fprintf(['Starting run ',num2str(run),' with main_run \n'])

%% Step 2: assign basics

% check that lambda is correctly given for important special case
checklambda(G,B,noise_cov,R,lambda);

% assign basics
[noise_cov,R] = makebasicspectral(noise_cov,R,G,B);

%% Step 3: calculate differentiation operators and finite difference matrices

fprintf('Creating differential operators ...\n');
[n,grid,region,D,D2,D4,fd1,fd2] = makediffopspectral4(bdim,n,bcon,region);


%% Step 4: calculate dynamics

fprintf('Creating dynamics operator ...\n');
%create ktensors
[fTens,GTens,BTens,noise_covTens,qTens,RTens] = maketensdyn(f,G,B,noise_cov,q,R,x,grid);

%create MATLAB functions
[fFunc,GFunc,BFunc,noise_covFunc,qFunc,RFunc] = makefuncdyn(f,G,B,noise_cov,q,R,x);

%% Step 5: calculate operator

fprintf('Creating PDE operator ...\n');
[op,conv,diff] = makeop(fTens,BTens,noise_covTens,qTens,D,D2,D4,lambda);

%% Step 6: add artificial diffusion to operator

if isempty(artdiff_option) == 0
    % add artificial diffusion
    h = [];
    [op] = addartdiff(op,conv,diff,artdiff_option,h,D2,fTens);  %%% TODO: NEED TO CHECK
end

%% Step 7: create boundary conditions

% create scaling for bc
if isempty(bsca) == 1 || isempty(regsca) == 1
    
    if sca_ver == 1
        [bscat, regscat] = make_bc_sca_var(op,grid,region,bcon);    
    elseif sca_ver == 2
        [bscat, regscat] = make_bc_sca(op,bcon,region,regval,als_options,fd1,grid,x,n);
    elseif sca_ver == 3
        % temporarily option made for comparing with old results
        op = (h(1)^2*h(2)^2)*op;
        bscat = ones(d,2);
        regscat = 1;
    else
        error('wrong specification on boundary scaling');
    end
    
end

if isempty(bsca)
    bsca = bscat;
end

if isempty(regsca)
    regsca = regscat;
end

% make bc
[bc] = makebc(bcon,bsca,grid,x,n);

%% Step 8: set up boundary conditions for operator
[op] = makebcopspectral(op,bcon,bsca,n,fd1);

%% Step 9: incooporate region in boundary conditions and operator

if isempty(region) == 0
    [op,bc] = incorpregion(op,bc,region,grid,regval,regsca);
end

%% Step 10: compress operator

op_uncomp = op; %save uncompressed op

fprintf('Attempt to compress operator, rank(op)=%d\n', ncomponents(op));
rank_op_uncomp = ncomponents(op);
tic;
start_compress = tic;

[op, err_op, iter_op, enrich_op, t_step_op, cond_op, noreduce] = als2(op,tol_err_op);

compress_time = toc(start_compress);
toc;
rank_op_comp = ncomponents(op);
fprintf('Number of components after compression, %d\n', ncomponents(op));

%% Step 11: solve system

disp('Beginning Solving');
tic;
start_solve = tic;

if isempty(als_variant) %original    
    [F, err, iter, Fcond, enrich, t_step, illcondmat, maxit, maxrank, F_cell, B_cell, b_cell] = ...
        als_sys(op,bc,[],tol_err,als_options,debugging);
    restart = []; %no restarts for original
    %save('F','F') %save just incase something does not work later
else %variant
    [F, err, iter, Fcond, enrich, t_step, illcondmat, maxit, maxrank, F_cell, B_cell, b_cell, restart] = ...
        als_sys_var(op,bc,[],tol_err,als_options,als_variant,debugging);
    %save('F','F') %save just incase something does not work later
end

time_solve = toc(start_solve);
toc;
disp('Solution complete');

%% Step 12: visualize results

if plotdata
% arrange input data
    F = arrange(F);
    plotsolve = {F,err,enrich,t_step,Fcond,grid};
    plotcomp = {op,err_op,enrich_op,t_step_op,cond_op};
    plotdebug = {F_cell,b_cell,B_cell};

% plot results from run als and compress operator
    try
        fprintf('Plotting results \n')
        visres(plotsolve,plotcomp,plotdebug,n,debugging,0,0,restart,run)
        fprintf('Plotting complete \n')
    catch
        fprintf('Could not visualize results \n')
    end
end

%% Step 13. run simulations

% % arrange input data
% sim_data = {lambda,grid,R,noise_cov,F,D,fFunc,GFunc,BFunc,qFunc,bdim,bcon,region};
% 
% % run simulation
% if isempty(sim_config) == 0
%     try
%         fprintf('Starting simulations \n')
%         sim_run(sim_config,sim_data,saveplots,savedata,run)
%         fprintf('Simulations complete \n')
%     catch
%         fprintf('Could not run simulation \n')
%     end
% end

%% Step 14: save data

time_whole = toc(start_whole);

if savedata == 1
    
    if debugging == 1
        
        % try to save F_cell,B_cell,b_cell separately since big files.
        try
            %save F_cell
            pathname = fileparts('./saved_data/');
            matfile = fullfile(pathname,['F_cell_run',num2str(run)]);
            save(matfile,'F_cell');
        catch
            fprintf('Could not save F_cell. Consider runs generating less data.');
        end
        try
            %save B_cell
            pathname = fileparts('./saved_data/');
            matfile = fullfile(pathname,['B_cell_run',num2str(run)]);
            save(matfile,'B_cell');
        catch
            fprintf('Could not save B_cell. Consider runs generating less data.');
        end
        try
            %save b_cell
            pathname = fileparts('./saved_data/');
            matfile = fullfile(pathname,['b_cell_run',num2str(run)]);
            save(matfile,'b_cell');
        catch
            fprintf('Could not save b_cell. Consider runs generating less data.');
        end
        
    end
    
    clear F_cell B_cell b_cell
    
    % save rest of the data
%     pathname = fileparts(save_folder);
%     matfile = fullfile(pathname,[save_folder,'spectral_rundata_run_',num2str(run)]);
    save([save_folder,'spectral_rundata_run_',num2str(run)]);
    
end

fprintf(['Run ',num2str(run),' with main_run is complete \n'])
