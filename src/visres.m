function visres(plotsolve,plotcomp,plotdebug,n,debugging,saveplots,savedata,restart,run)
% VISRES plots results from a run by MAIN_RUN.
% Inputs:
%   plotsolve - cell with the following data:
%      F - obtained solution from als run.
%      err - achieved error in als run.
%      enrich - when F increased in rank during als run.
%      t_step - how much time iterations in als run took.
%      Fcond - contition number of F during als run.
%      grid - grid points in each dimension.
%   plotcompress - cell with the following data:
%      op - compressed operator
%      err_op - achieved error during compress operator.
%      enrich_op - when op increased in rank during compress operator.
%      t_step_op - how much time iterations took when compressing operator.
%      cond_op - condition number of op during compress operator
%   plotdebug is a cell with the following data:
%      F_cell - F during als run.
%      b_cell - RHS vectors during als run.
%      B_cell - als matrices during als run.
%   n - number of grid in each dimension
%      debugging - 1 if F_cell,B_cell,b_cell shall be documented.
%      saveplots - 1 if plots shall be saved.
%      run - index separating different runs from MAIN_RUN
% Outputs: (OBS: are just saved)
%   visres_output - cell with the following data:
%      (np = no precondition, p = precondition)
%      normF_vector_np - norms of F during als run with np.
%      normF_vector_p  - norms of F during als run with p.
%      normb_vector_np - norms of RHS vector during als run with np.
%      normb_vector_p  - norms of RHS vector during als run with p.
%      condB_vector_np - condition number for als matrices during als run 
%      with np.
%      condB_vector_p  - condition number for als matrices during als run
%      with p.
%      condB_each_dim  - condition number for als matrices during als run. 
%      Row is iteration number and column is dimension.
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 2015

%% Extract data
[F,err,enrich,t_step,Fcond,grid] = deal(plotsolve{:});

[op,err_op,enrich_op,t_step_op,cond_op] = deal(plotcomp{:});

[F_cell,b_cell,B_cell] = deal(plotdebug{:});

d = length(grid);

%% Plot data from als run

% general data
fig = figure;

% convergence
subplot(2,2,1);
semilogy(err,'LineWidth',2);
xlabel('iteration','FontSize',14);
ylabel('error','FontSize',14);
title('Convergence of F','FontSize',14);
set(gca,'FontSize',12);
hold on;
for i=1:length(enrich)
    plot(enrich(i),err(enrich(i)),'r*','MarkerSize',10);
end
for i=1:length(restart)
    plot(restart(i),err(restart(i)),'gs','MarkerSize',10);
end
hold off;

% time per iteration
subplot(2,2,2);
plot(t_step,'LineWidth',2);
xlabel('iteration','FontSize',14);
ylabel('Time(s)','FontSize',14);
title('Computation Time','FontSize',14);
set(gca,'FontSize',12)
hold on;
for i=1:length(enrich)
    plot(enrich(i),t_step(enrich(i)),'r*','MarkerSize',10);
end
for i=1:length(restart)
    plot(restart(i),err(restart(i)),'gs','MarkerSize',10);
end
hold off;

% weighting of basis functions
subplot(2,2,3);
semilogy(F.lambda,'LineWidth',2);
xlabel('Basis Number','FontSize',14);
ylabel('Weight','FontSize',14);
title('Composition of F','FontSize',14);
set(gca,'FontSize',12)

% conditioning number for F
subplot(2,2,4);
plot(Fcond,'LineWidth',2);
xlabel('iteration','FontSize',14);
ylabel('Cond','FontSize',14);
title('Condition number of F','FontSize',14);
set(gca,'FontSize',12)

if saveplots == 1
    pathname = fileparts('./saved_data/');
    matfig = fullfile(pathname,['solve_run',num2str(run)]);
    saveas(fig,matfig);
end

% additional plots if dimension equals 1 or 2
if d == 1
    
    sol = double(F)';
    
    fig = figure;
    plot(grid{1},sol)
    xlabel('x'); ylabel('\Psi')
    title('\Psi');
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['solution_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
end

if d == 2
    
    sol = double(F)';
    [x1mesh,x2mesh] = meshgrid(grid{1},grid{2});
    
    figure;
    contour(grid{2},grid{1},reshape(sol,n(1),n(2)))
    xlabel('x_1');ylabel('x_2');
    
    % plot solution as 2D-plot
    fig = figure;
    imagesc(grid{1},grid{2},sol); axis xy
    xlabel('x_1');ylabel('x_2');
    title('\Psi - 2D');
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['solution_2D_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
    % plot solution as 3D-plot
    fig = figure;
    surf(x1mesh,x2mesh,sol,'EdgeColor','none');
    xlabel('x_1');ylabel('x_2');zlabel('\Psi');
    title('\Psi - 3D');
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['solution_3D_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
end

%% Plot data from compressing operator

fig = figure;

% convergence
subplot(2,2,1);
semilogy(err_op,'LineWidth',2);
xlabel('iteration','FontSize',14);
ylabel('error','FontSize',14);
title('Convergence of Operator','FontSize',14);
set(gca,'FontSize',12);
hold on;
for i=1:length(enrich_op)
    plot(enrich_op(i),err_op(enrich_op(i)),'r*','MarkerSize',10);
end
hold off;

% time per iteration
subplot(2,2,2);
plot(t_step_op,'LineWidth',2);
xlabel('iteration','FontSize',14);
ylabel('Time(s)','FontSize',14);
title('Computation Time','FontSize',14);
set(gca,'FontSize',12)
hold on;
for i=1:length(enrich_op)
    plot(enrich_op(i),t_step_op(enrich_op(i)),'r*','MarkerSize',10);
end
hold off;

% weighting of basis functions
subplot(2,2,3);
semilogy(op.lambda,'LineWidth',2);
xlabel('Basis Number','FontSize',14);
ylabel('Weight','FontSize',14);
title('Composition of Operator','FontSize',14);
set(gca,'FontSize',12)

if saveplots == 1
    pathname = fileparts('./saved_data/');
    matfig = fullfile(pathname,['compress_op_run',num2str(run)]);
    saveas(fig,matfig);
end

% conditioning number for op
subplot(2,2,4);
plot(cond_op,'LineWidth',2);
xlabel('iteration','FontSize',14);
ylabel('Cond','FontSize',14);
title('Condition number of Operator','FontSize',14);
set(gca,'FontSize',12)

if saveplots == 1
    pathname = fileparts('./saved_data/');
    matfig = fullfile(pathname,['compress_op_run',num2str(run)]);
    saveas(fig,matfig);
end

%% plot first five basis functions in each dimension
pdim = ceil(sqrt(d));
fig = figure;

if length(enrich) < 5
    
    for i=1:d
        
        subplot(pdim,pdim,i)
        tit = sprintf('dimension = %d', i);
        %fil = sprintf('basis_%d.eps', i);
        plot(F.U{i}(:,1:length(enrich)),'LineWidth',2);
        hold on;
        plot([1 n(i)],[0 0], '--', 'Color',[.3 .3 .3], 'LineWidth',1);
        hold off;
        title(tit,'FontSize',12);
        axis([1 n(i) -inf inf])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        axis off
        
    end
    
else
    
    for i=1:d
        
        subplot(pdim,pdim,i)
        tit = sprintf('dimension = %d', i);
        %fil = sprintf('basis_%d.eps', i);
        plot(F.U{i}(:,1:5),'LineWidth',2);
        hold on;
        plot([1 n(i)],[0 0], '--', 'Color',[.3 .3 .3], 'LineWidth',1);
        hold off;
        title(tit,'FontSize',12);
        axis([1 n(i) -inf inf])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        axis off;
        
    end
    
end

if saveplots == 1
    pathname = fileparts('./saved_data/');
    matfig = fullfile(pathname,['basis_first_run',num2str(run)]);
    saveas(fig,matfig);
end

%% plot all basis functions in each dimension
fig = figure;

for i=1:d
    
    subplot(pdim,pdim,i)
    tit = sprintf('dimension = %d', i);
    %fil = sprintf('basis_%d.eps', i);
    plot(F.U{i}(:,1:end),'LineWidth',2);
    hold on;
    plot([1 n(i)],[0 0], '--', 'Color',[.3 .3 .3], 'LineWidth',1);
    hold off;
    title(tit,'FontSize',12);
    axis([1 n(i) -inf inf])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis off;
    
end

if saveplots == 1
    pathname = fileparts('./saved_data/');
    matfig = fullfile(pathname,['basis_all_run',num2str(run)]);
    saveas(fig,matfig);
end

%% plot norm of solutions F
if debugging == 1
    
    sizeFcell = size(F_cell);
    
    % no precondition
    Fcell_vector_np = F_cell{1,4};
    Fcell_iter_np = [];
    Fcell_enrich_np = 1;
    Fcell_restart_np = [];
    
    % precondition
    Fcell_vector_p = F_cell{1,4};
    Fcell_iter_p = [];
    Fcell_enrich_p = 1;
    
    % create Fcell_vector from Fcell
    for i = 1:sizeFcell(1)
        
        % OBS: adding rank and new iteration do not coincide for F
        Fcell_iter_np = [Fcell_iter_np, length(Fcell_vector_np)+1];
        Fcell_iter_p = [Fcell_iter_p, length(Fcell_vector_p)+1];
        
        c = F_cell{i,1};
        Fcell_vector_np = [Fcell_vector_np, c];
        Fcell_vector_p = [Fcell_vector_p, c];
        
        if any(restart==i)
            Fcell_restart_np = [Fcell_restart_np length(Fcell_vector_np)+1];
        end
        
        c = F_cell{i,2};
        if isempty(c) == 0
            Fcell_vector_p = [Fcell_vector_p, c];
        end
        
        c = F_cell{i,3};
        if isempty(c) == 0
            Fcell_vector_np = [Fcell_vector_np, c];
            Fcell_vector_p = [Fcell_vector_p, c];
            
            Fcell_enrich_np = [Fcell_enrich_np, length(Fcell_vector_np)];
            Fcell_enrich_p = [Fcell_enrich_p, length(Fcell_vector_p)];
        end
        
    end
    
    % take norm of F_cell_vector
    normF_vector_np = ones(1,length(Fcell_vector_np));
    normF_vector_p = ones(1,length(Fcell_vector_p));
    
    for i = 1:length(Fcell_vector_np)
        normF_vector_np(i) = norm(Fcell_vector_np{i});
    end
    
    for i = 1:length(Fcell_vector_p)
        normF_vector_p(i) = norm(Fcell_vector_p{i});
    end
    
    % plot without precondition
    fig = figure;
    semilogy(normF_vector_np,'LineWidth',1);
    xlabel('updated F','FontSize',14);
    ylabel('norm(F)','FontSize',14);
    title('norm(F)','FontSize',14);
    set(gca,'FontSize',12)
    hold on
    for i = 1:length(Fcell_iter_np)
        plot(Fcell_iter_np(i),normF_vector_np(Fcell_iter_np(i)),'bo','MarkerSize',5);
    end
    for i = 1:length(Fcell_enrich_np)
        plot(Fcell_enrich_np(i),normF_vector_np(Fcell_enrich_np(i)),'r*','MarkerSize',5);
    end
    for i = 1:length(Fcell_restart_np)
        plot(Fcell_restart_np(i),normF_vector_np(Fcell_restart_np(i)),'gs','MarkerSize',5);
    end
    hold off;
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['normF_np_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
    % plot with precondition
    fig = figure;
    semilogy(normF_vector_p,'LineWidth',1);
    xlabel('updated F','FontSize',14);
    ylabel('norm(F)','FontSize',14);
    title('norm(F) - including precondition','FontSize',14);
    set(gca,'FontSize',12)
    hold on
    for i = 1:length(Fcell_iter_p)
        plot(Fcell_iter_p(i),normF_vector_p(Fcell_iter_p(i)),'bo','MarkerSize',5);
    end
    for i = 1:length(Fcell_enrich_p)
        plot(Fcell_enrich_p(i),normF_vector_p(Fcell_enrich_p(i)),'r*','MarkerSize',5);
    end
    hold off;
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['normF_p_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
end

%% plot norm of b

if debugging == 1
    
    size_bcell = size(b_cell);
    
    % no precondition
    normb_vector_np = [];
    iter_vector_np = 1;
    enrich_vector_np = 1;
    restart_vector_np = [];
    
    % precondition
    normb_vector_p = [];
    iter_vector_p = 1;
    enrich_vector_p = 1;
    
    for i = 1:size_bcell(1)
        
        c = b_cell{i,1};
        
        iter_vector_np = [iter_vector_np length(normb_vector_np)+1];
        iter_vector_p = [iter_vector_p length(normb_vector_p)+1];
        
        add = ones(1,length(c));
        for j = 1:length(c)
            add(j) = norm(c{j});
        end
        normb_vector_np = [normb_vector_np add];
        normb_vector_p = [normb_vector_p add];
        
        if any(restart==i)
            restart_vector_np = [restart_vector_np length(normb_vector_np)+1];
        end
        
        if isempty(b_cell{i,2}) == 0 && i ~= size_bcell(1)
            enrich_vector_np = [enrich_vector_np length(normb_vector_np)+1];
            enrich_vector_p = [enrich_vector_p length(normb_vector_p)+1];
            
            c = b_cell{i,2};
            
            add = ones(1,length(c));
            for j = 1:length(c)
                add(j) = norm(c{j});
            end
            normb_vector_p = [normb_vector_p add];
            
        end
    end
    
    % plot without precondition
    fig = figure;
    semilogy(normb_vector_np,'LineWidth',1);
    xlabel('updated RHS vector','FontSize',14);
    ylabel('norm','FontSize',14);
    title('norm for RHS vectors','FontSize',14);
    set(gca,'FontSize',12)
    hold on
    % add iter
    for i = 1:length(iter_vector_np)
        plot(iter_vector_np(i), normb_vector_np(iter_vector_np(i)),'bo','MarkerSize',5);
    end
    % add enrich
    for i = 1:length(enrich_vector_np)
        plot(enrich_vector_np(i), normb_vector_np(enrich_vector_np(i)),'r*','MarkerSize',5);
    end
    % add restart
    for i = 1:length(restart_vector_np)
        plot(restart_vector_np(i), normb_vector_np(restart_vector_np(i)),'gs','MarkerSize',5);
    end
    hold off
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['normb_np_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
    % plot with precondition
    fig = figure;
    semilogy(normb_vector_p,'LineWidth',1);
    xlabel('updated RHS vector','FontSize',14);
    ylabel('norm','FontSize',14);
    title('Norm for RHS vectors - including precondition','FontSize',14);
    set(gca,'FontSize',12)
    hold on
    % add iter
    for i = 1:length(iter_vector_p)
        plot(iter_vector_p(i), normb_vector_p(iter_vector_p(i)),'bo','MarkerSize',5);
    end
    % add enrich
    for i = 1:length(enrich_vector_p)
        plot(enrich_vector_p(i), normb_vector_p(enrich_vector_p(i)),'r*','MarkerSize',5);
    end
    hold off
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['normb_p_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
end

%% plot square root of norm of ALS matrix without precondition.
% By taking square root, you are able to compare the magnitude with the
% magnitude of the RHS vector b

if debugging == 1
    
    size_Bcell = size(B_cell);
    
    % no precondition
    normB_vector_np = [];
    iter_vector_np = 1;
    enrich_vector_np = 1;
    restart_vector_np = [];
    
    for i = 1:size_bcell(1)
        
        c = B_cell{i,1};
        
        iter_vector_np = [iter_vector_np length(normB_vector_np)+1];
        
        add = ones(1,length(c));
        for j = 1:length(c)
            add(j) = sqrt(norm(c{j}));
        end
        normB_vector_np = [normB_vector_np add];
        
        if any(restart==i)
            restart_vector_np = [restart_vector_np length(normB_vector_np)+1];
        end
        
        if isempty(B_cell{i,2}) == 0 && i ~= size_bcell(1)
            enrich_vector_np = [enrich_vector_np length(normB_vector_np)+1];
        end
    end
    
    % plot
    fig = figure;
    semilogy(normB_vector_np,'LineWidth',1);
    xlabel('updated ALS matrix','FontSize',14);
    ylabel('square root of norm','FontSize',14);
    title('Square root of the ALS matrix norm','FontSize',14);
    set(gca,'FontSize',12)
    hold on
    % add iter
    for i = 1:length(iter_vector_np)
        plot(iter_vector_np(i), normB_vector_np(iter_vector_np(i)),'bo','MarkerSize',5);
    end
    % add enrich
    for i = 1:length(enrich_vector_np)
        plot(enrich_vector_np(i), normB_vector_np(enrich_vector_np(i)),'r*','MarkerSize',5);
    end
    % add restart
    for i = 1:length(restart_vector_np)
        plot(restart_vector_np(i), normB_vector_np(restart_vector_np(i)),'gs','MarkerSize',5);
    end
    hold off
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['normB_np_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
end

%% plot condition number of ALS matrices with/without precondition and in each dimension
if debugging == 1
    
    size_Bcell = size(B_cell);
    
    % no precondition
    condB_vector_np = [];
    iter_vector_np = 1;
    enrich_vector_np = 1;
    restart_vector_np = [];
    
    % precondition
    condB_vector_p = [];
    iter_vector_p = 1;
    enrich_vector_p = 1;
    
    % each dimension plots
    enrich_each_dim = [];
    restart_each_dim = [];
    
    for i = 1:size_Bcell(1)
        
        c = B_cell{i,1};
        
        iter_vector_np = [iter_vector_np length(condB_vector_np)+1];
        iter_vector_p = [iter_vector_p length(condB_vector_p)+1];
        
        add = ones(1,length(c));
        for j = 1:length(c)
            add(j) = cond(c{j});
        end
        condB_vector_np = [condB_vector_np add];
        condB_vector_p = [condB_vector_p add];
        
        if any(restart==i)
            restart_vector_np = [restart_vector_np length(condB_vector_np)+1];
            restart_each_dim = [restart_each_dim i];
        end
        
        if isempty(B_cell{i,2}) == 0 && i ~= size_Bcell(1)
            enrich_vector_np = [enrich_vector_np length(condB_vector_np)+1];
            enrich_vector_p = [enrich_vector_p length(condB_vector_p)+1];
            enrich_each_dim = [enrich_each_dim i];
            
            c = B_cell{i,2};
            
            add = ones(1,length(c));
            for j = 1:length(c)
                add(j) = cond(c{j});
            end
            condB_vector_p = [condB_vector_p add];
            
        end
    end
    
    condB_each_dim = reshape(condB_vector_np,length(condB_vector_np)/d,d);
    
    % plot without precondition
    fig = figure;
    semilogy(condB_vector_np,'LineWidth',1);
    xlabel('updated ALS matrix','FontSize',14);
    ylabel('condition number','FontSize',14);
    title('Condition number for ALS matrices','FontSize',14);
    set(gca,'FontSize',12)
    hold on
    % add iter
    for i = 1:length(iter_vector_np)
        plot(iter_vector_np(i), condB_vector_np(iter_vector_np(i)),'bo','MarkerSize',5);
    end
    % add enrich
    for i = 1:length(enrich_vector_np)
        plot(enrich_vector_np(i), condB_vector_np(enrich_vector_np(i)),'r*','MarkerSize',5);
    end
    % add restart
    for i = 1:length(restart_vector_np)
        plot(restart_vector_np(i), condB_vector_np(restart_vector_np(i)),'gs','MarkerSize',5);
    end
    hold off
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['condB_np_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
    % plot with precondition
    fig = figure;
    semilogy(condB_vector_p,'LineWidth',1);
    xlabel('updated ALS matrix','FontSize',14);
    ylabel('condition number','FontSize',14);
    title('Condition number for ALS matrices - including precondition','FontSize',14);
    set(gca,'FontSize',12)
    hold on
    % add iter
    for i = 1:length(iter_vector_p)
        plot(iter_vector_p(i), condB_vector_p(iter_vector_p(i)),'bo','MarkerSize',5);
    end
    % add enrich
    for i = 1:length(enrich_vector_p)
        plot(enrich_vector_p(i), condB_vector_p(enrich_vector_p(i)),'r*','MarkerSize',5);
    end
    hold off
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['condB_p_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
    % plot in each dimension (without precondition)
    pdim = ceil(sqrt(d));
    fig = figure;
    
    for i = 1:d
        subplot(pdim,pdim,i);
        condB_each_dim_i = condB_each_dim(:,i);
        semilogy(condB_each_dim_i,'LineWidth',1);
        xlabel('updated ALS matrix','FontSize',8);
        ylabel('condition number','FontSize',8);
        title(['Condition number for ALS matrices in dimension ',num2str(i)],'FontSize',8);
        set(gca,'FontSize',12)
        hold on
        % add enrich
        for j = 1:length(enrich_each_dim)
            plot(enrich_each_dim(j), condB_each_dim_i(enrich_each_dim(j)),'r*','MarkerSize',5);
        end
        % add restart
        for j = 1:length(restart_each_dim)
            plot(restart_each_dim(j), condB_each_dim_i(restart_each_dim(j)),'gs','MarkerSize',5);
        end
        hold off
    end
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['condB_each_dim_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
end

%% save data
if savedata == 1
    clear fig;
    pathname = fileparts('./saved_data/');
    matfile = fullfile(pathname,['visdata_run',num2str(run)]);
    save(matfile);
end

