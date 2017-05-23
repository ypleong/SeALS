function [fTens] = fcell2ftens(fcell,grid)
%FUNC2TENS Creates a function in ktensor form from cell array
%representation.
% Input:
%   fcell - the function in cell array representation.
%   grid - the discretization grid.
% Output:
%   fTens - the function as a cell array. fTens{i,j} is the ktensor of
%   element (i,j) in the original function f.
%
% See also FSYM2FCELL, MAKETENSDYN.

% Elis Stefansson, Aug 5 2015
% The function will be denoted f below.

size_fcell = size(fcell);
fTens = cell(size_fcell);

for i = 1:size_fcell(1)
    
    for j = 1:size_fcell(2);
        
        c = fcell{i,j}; % cell repr. of the function in entry (i,j) of f.
        
        for ii = 1:length(c)
            
            cc = c{ii};
            size_cc = size(cc);
            
            % vectors for the corresponding ktensor
            U = cell(1,length(grid));
            
            for jj = 1:size_cc(1)
                
                % take function and corresponding dimension
                c_func = cc{jj,1};
                dim = cc{jj,2};
                c_grid = grid{dim};
                c_vect = c_func(c_grid); %convert to vector
                
                % if c_func is constant, length(c_vect) = 1
                lgd = length(grid{dim});
                if length(c_vect) ~= lgd
                    c_vect = c_vect*ones(lgd,1);
                end
                
                if isempty(U{dim}) == 1
                    U{dim} = c_vect;
                else
                    U{dim} = U{dim}.*c_vect; %multiply terms in same dim
                end
                
            end
            
            % add ones to dimension not included in U already
            for k = 1:length(grid)
                if isempty(U{k}) == 1
                    U{k} = ones(length(grid{k}),1);
                end
            end
            
            if isempty(fTens{i,j})
                fTens{i,j} = fixsigns(arrange(ktensor(U)));
            else
                fTens{i,j} = fixsigns(arrange(fTens{i,j}+ktensor(U)));
            end
            
        end
    end
end

end