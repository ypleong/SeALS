function sim_run(sim_config,sim_data,saveplots,savedata,run)
% SIMULATION simulates trajectories starting from x0_list using results
% from MAIN_RUN.
% Inputs:
%   sim_config is cell with the following data:
%      T - end time for simulations.
%      dt - time step.
%      x0_list - matrix with starting points; each row corresponds to a
%      starting point.
%      new_noise_cov - variance vector for noise in simulation. The noise
%      from main_run_data will be used if new_noise_cov is empty.
%      new_region - goal region in simulation, replaces the old region. The
%      old region will be used if new_region is empty.
%      run - run index.
%   sim_data is a cell with the following data:
%      lambda - the scalar lambda.
%      grid - the discretization grid.
%      R - the cost matrix.
%      noise_cov - the covariance matrix of the noise.
%      F - the solution F (desirability function) as ktensor.
%      D - cell with first order differentiation operators. D{i} is
%      differentiation with respect to dimension i.
%      fFunc - f as MATLAB function.
%      GFunc - G as MATLAB function.
%      BFunc - B as MATLAB function.
%      qFunc - q as MATLAB function.
%      bdim - the dimensions of the hyperrectangle domain.
%      bcon - the boundary conditions. bcon{i} is boundary conditions in
%      dimension i:
%         bcon{i} = {'p'}. Periodic
%         bcon{i} = {'d',val_lo,val_up}. Dirichlet with val_lo and val_up 
%         for lower and upper boundary in dimension i.
%         bcon{i} = {'n',val_lo,val_up}. Neumann with val_lo and val_up for
%         lower and upper boundary in dimension i.
%         where all the values are for the desirability function.
%      region - the dimensions of the goal region.
%   saveplots - 1 if plots will be saved.
%   savedata - 1 if data will be saved.
%   run - run index.
%   See also MAIN_PROGRAM.
% Outputs: (OBS: only saved!)
%   simulation plots - state trajectories, control trajectories and cost
%   trajectories.
%   simulation data - saved as sim_data_run where run is the corresponding
%   run index.
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 2015

%% extract data
[T,dt,x0_list,new_noise_cov,new_region] = deal(sim_config{:});

[lambda,grid,R,noise_cov,F,D,fFunc,GFunc,BFunc,qFunc,bdim,bcon,region] = deal(sim_data{:});

d = length(grid);
ninputs = length(R);

if isempty(new_noise_cov) == 0
    sigma = new_noise_cov; %use other noise
else
    sigma = diag(noise_cov); %convert cov-matrix to vector
end

if isempty(new_region) == 0
    region = new_region; %use other goal region
end

dF = cell(d,1);
for i = 1:d
    dF{i} = SRMultV(D{i},F);
end

%% simulation

size_x0_list = size(x0_list);
n_trial = size_x0_list(2);
t = 0:dt:T;
t_length = length(t);

traj = zeros(n_trial,t_length,d);
traj_u = zeros(n_trial,t_length,ninputs);
c_traj = zeros(n_trial,t_length);

end_list = t_length*ones(n_trial);
periodic_list = cell(n_trial,d); %index for periodic jumps for each dim
conv_list = zeros(n_trial,1);
end_point = zeros(n_trial,d);

fprintf('Simulation number = \n');

for m = 1:n_trial
    
    current = x0_list(:,m);
    traj(m,1,:) = current;
    
    periodic_list_m = cell(1,d);
    for i = 1:d
        periodic_list_m{i} = 0;
    end
    
    for k = 1:t_length-1
        
        % calculate current values
        c_f = fFunc(current);
        c_G = GFunc(current);
        c_B = BFunc(current);
        c_F = EvalT(F,current,grid);
        c_dF = EvalTMat(dF,current,grid);
        c_q = qFunc(current);
        
        c_n = normrnd(0,sigma);
        c_u = (1/c_F)*lambda*inv(R)*c_G'*c_dF;
        
        next = current + dt*(c_f+c_G*c_u+c_B*c_n);
        cost = c_q + 0.5*c_u'*R*c_u;
        current = next;
        
        traj_u(m,k,:) = c_u;
        traj(m,k+1,:) = current;
        c_traj(m,k) = cost;
        
        % check where current state is
        region_check = zeros(1,d);
        
        for i=1:d
            
            bdim_i = bdim(i,:);
            bdim_i_l = bdim_i(2)-bdim_i(1);
            bcon_i = bcon{i};
            
            if current(i) < bdim_i(1)
                
                if length(bcon_i) == 1 %then periodic
                    current(i) = current(i)+bdim_i_l;
                    periodic_list_m{i} = [periodic_list_m{i} k+1];
                else
                    end_list(m) = k+1;
                    conv_list(m) = -1; %goes unstable
                    for j = 1:d
                        periodic_list_m{j} = [periodic_list_m{j} k+1];
                    end
                end
                
            elseif current(i) > bdim_i(2)
                
                if length(bcon_i) == 1 %then periodic
                    current(i) = current(i)-bdim_i_l;
                    periodic_list_m{i} = [periodic_list_m{i} k+1];
                else
                    end_list(m) = k+1;
                    conv_list(m) = -1;
                    for j = 1:d
                        periodic_list_m{j} = [periodic_list_m{j} k+1];
                    end
                end
                
            end
            
            if isempty(region) == 0
                if current(i) > region(i,1) && current(i) < region(i,2)
                    region_check(i) = 1;
                end
            end
            
        end
        
        if conv_list(m) == -1
            break;
        end
        
        if isempty(region) == 0
            if sum(region_check) == d
                end_list(m) = k+1;
                conv_list(m) = 1; %stable to goal region
                for j=1:d
                    periodic_list_m{j} = [periodic_list_m{j} k+1];
                end
                break;
            end
        end
        
    end
    
    end_point(m,:) = current;
    for j = 1:d
        periodic_list_m{j} = [periodic_list_m k+1];
        periodic_list_m{j} = unique(cell2mat(periodic_list_m{j}));
        periodic_list{m,j} = periodic_list_m{j};
    end
    fprintf(' %d',m);
    
end

fprintf('\n');

%% calculate total cost
c_tot = dt*sum(c_traj,2);

%% plot state trajectories

pdim = ceil(sqrt(d));
ccc = jet(n_trial);

% plot state trajectories for each dimension seperately
fig = figure;
for i=1:d
    
    subplot(pdim,pdim,i);
    hold on
    for j = 1:n_trial
        
        c_per_list = periodic_list{j,i};
        
        for k = 1:length(c_per_list)-1
            
            cc_per_list = c_per_list(k)+1:c_per_list(k+1);
            plot( t(cc_per_list), squeeze( traj(j,cc_per_list,i) ), 'color', ccc(j,:) );
            
        end
        
    end
    
    if isempty(region) == 0
        plot([t(1) t(end_list(j))],[region(i,1) region(i,1)], '--', 'Color',[.3 .3 .3], 'LineWidth',1);
        plot([t(1) t(end_list(j))],[region(i,2) region(i,2)], '--', 'Color',[.3 .3 .3], 'LineWidth',1);
    end
    hold off
    xlabel('time(s)')
    axis([0 inf bdim(i,1) bdim(i,2)])
    title(['state traj for x_',num2str(i)])
    
end

if saveplots == 1
    pathname = fileparts('./saved_data/');
    matfig = fullfile(pathname,['traj_plot_run',num2str(run)]);
    saveas(fig,matfig);
end

% special case when d = 2, trajectories in one plot.
if d == 2
    
    % combine peridic_list for dimension 1 and 2
    per_list_d2 = cell(1,n_trial);
    for j = 1:n_trial
        per_list_d2{j} = unique([periodic_list{j,1} periodic_list{j,2}]);
    end
    
    fig = figure;
    hold on
    for j = 1:n_trial
        
        c_per_list = per_list_d2{j};
        
        for k = 1:length(c_per_list)-1
            
            cc_per_list = c_per_list(k)+1:c_per_list(k+1);
            plot(squeeze((traj(j,cc_per_list,1))),squeeze((traj(j,cc_per_list,2))),'color',ccc(j,:))
            
        end
        
    end
    
    if isempty(region) == 0
        rectangle('Position',[region(1,1),region(2,1),region(1,2)-region(1,1),region(2,2)-region(2,1)],...
            'Curvature',[0,0],'FaceColor','r')
    end
    xlabel('x_1')
    ylabel('x_2')
    title('state trajectories')
    axis([bdim(1,1) bdim(1,2) bdim(2,1) bdim(2,2)])
    set(gca,'fontsize', 16)
    box on
    hold off;
    
    if saveplots == 1
        pathname = fileparts('./saved_data/');
        matfig = fullfile(pathname,['traj_2Dplot_run',num2str(run)]);
        saveas(fig,matfig);
    end
    
end

%% plot control trajectories

pdim = ceil(sqrt(ninputs));
ccc = jet(n_trial);

% plot control trajectories for each control dimension seperately
fig = figure;
for i=1:ninputs
    
    subplot(pdim,pdim,i);
    hold on
    for j = 1:n_trial
        plot( t(1:end_list(j)-1), squeeze( traj_u(j,1:end_list(j)-1,i) ), 'color', ccc(j,:) );
    end
    hold off
    xlabel('time(s)')
    title(['control traj:s for u_',num2str(i)])
    
end

if saveplots == 1
    pathname = fileparts('./saved_data/');
    matfig = fullfile(pathname,['traj_u_plot_run',num2str(run)]);
    saveas(fig,matfig);
end

%% plot cost trajectory

ccc = jet(n_trial);
fig = figure;

hold on
for j = 1:n_trial
    plot( t(1:end_list(j)-1), c_traj(j,1:end_list(j)-1), 'color', ccc(j,:) );
end
hold off
xlabel('time(s)')
title('cost trajectories')

if saveplots == 1
    pathname = fileparts('./saved_data/');
    matfig = fullfile(pathname,['c_traj_plot_run',num2str(run)]);
    saveas(fig,matfig);
end

%% save data

if savedata == 1
    clear fig;
    pathname = fileparts('./saved_data/');
    matfile = fullfile(pathname,['simdata_run',num2str(run)]);
    save(matfile);
end



