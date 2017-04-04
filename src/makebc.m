function [bc] = makebc(bcon,bsca,grid,x,n)
% MAKEBC creates the boundary conditions with scaling.
% Input:
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
%   grid - the discreatization grid.
%   x - the symbolic state space vector.
%   n - the number of grid points in each dimension.
% Output:
%   bc - the boundary conditions as a ktensor
%
% Remark: The boundary points common for more than one dimension will 
% attain the lowest dimension boundary value.
%
% See also MAKE_BC_SCA, MAIN_RUN.

% Elis Stefansson, Aug 6 2015

d = length(n);
bc = [];

for i = 1:d %add boundary conditions to each dimension
    
    bconi = bcon{i};
    
    if strcmp(bconi{1},'n') || strcmp(bconi{1},'d') %not periodic
        
        
        % 1. create lower boundary tensor
        bcl = bconi{2};
        if isa(bcl,'sym') == 0 %not sym (constant)
            bcl = sym(bcl);
        end
        bcl_cell = fsym2fcell(bcl,x);
        bcl_tens = fcell2ftens(bcl_cell,grid);
        bcl_tens = bcl_tens{1}; %from cell to tensor
        
        % do not intefer with bc for lower dimensions
        for k = 1:(i-1)
            bcl_tens.U{k}(1,:) = 0;
            bcl_tens.U{k}(end,:) = 0;
        end
        
        % just the lower boundary points
        Ui = bcl_tens.U{i};
        si = size(Ui);
        newUi = zeros(si);
        newUi(1,:) = Ui(1,:);
        bcl_tens.U{i} = newUi;
        
        % scale the boundary
        bcl_tens = bsca(i,1)*bcl_tens;
        
        
        % 2. create upper boundary tensor
        bcu = bconi{3};
        if isa(bcu,'sym') == 0 %not sym (constant)
            bcu = sym(bcu);
        end
        bcu_cell = fsym2fcell(bcu,x);
        bcu_tens = fcell2ftens(bcu_cell,grid);
        bcu_tens = bcu_tens{1};
        
        % do not intefer with bc for lower dimensions
        for k = 1:(i-1)
            bcu_tens.U{k}(1,:) = 0;
            bcu_tens.U{k}(end,:) = 0;
        end
        
        % just the upper boundary points
        Ui = bcu_tens.U{i};
        si = size(Ui);
        newUi = zeros(si);
        newUi(end,:) = Ui(end,:);
        bcu_tens.U{i} = newUi;
        
        % scale the boundary
        bcu_tens = bsca(i,1)*bcu_tens;
        
        % 3. combine
        bci = bcl_tens+bcu_tens;
        
        if isempty(bc) == 1
            bc = bci;
        else
            bc = bc+bci;
        end
        
    end
    
end

if isempty(bc) == 1 %all periodic or vanishing
    
    bcU = cell(1,d);
    for i=1:d
        bcU{i} = zeros(n(i),1);
    end
    
    bc = ktensor(bcU);
    
end

