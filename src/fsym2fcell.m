function [fcell] = fsym2fcell(fsym,x)
% FSYM2FCELL Creates a function in cell array representation from symbolic
% representation.
% Inputs:
%   fsym - function in symbolic representation. Possibly a matrix function.
%   x - the symbolic vector.
% Outputs:
%   fcell - the function f as a cell array.
%
% Explanation of the cell array representation:
%
%   fcell{i,j} is the cell array representation of the (scalar) function in
%   entry (i,j) in fsym. For simplicity, denote the scalar symbolic 
%   function by f_sca and the cell representation by f_sca_cell. 
%   f_sca_cell{k} corresponds to a additive term in f_sca where
%   f_sca_cell{k}{l} corresponds to a multiplicatve term of just one 
%   variable in the additive term of f_sca.
%
%   Thus, the cell representation first splits
%   the components of the fsym, then splits up the
%   additive terms and then splits up the one-dimensional multiplicative 
%   terms.
%
%   Note: it is important that the multiplicative in fsym are 
%   one-dimensional since the conversion would not work otherwise. A not 
%   alowed example is fsym = (x(1)+x(2))*x(1). Consider 
%   fsym = x(1)*x(1)+x(1)*x(1) instead.
%
%   An example:
%   x = sym('x',[2,1]) %two dimensions
%   f = x(1)*x(2)+x(2)
%   becomes in cell array represenation
%   fcell = { { @(x2)x2, [2] }, { @(x1)x1, [1] ; @(x2)x2, [2]} }
%                  ^ matlab function corresponding to the one-dim mult term
%                         ^ and which dimension.
%
%
% See also MAKETENSDYN.

% Elis Stefansson, Aug 5 2015

size_fsym = size(fsym);
fcell = cell(size_fsym);

for i=1:size_fsym(1)
    
    for j=1:size_fsym(2)
        
        c_sym = fsym(i,j);
        ch = children(c_sym);
        
        c_check = ch(1);
        for k=2:length(ch)
            c_check = c_check+ch(k);
        end
        
        % check if we have addition terms
        if isequaln(c_check,c_sym)
            % case: several addition terms, x(1)+x(1)*x(2)^2 for
            % example, or just one term on the form x(i).
            
            % create cell for add. terms
            add_cell = cell(1,length(ch));
            
            for ii = 1:length(ch)
                
                if length(symvar(ch(ii))) < 2
                    % case: just one variable term or constant term
                    
                    % create cell for mult. terms
                    mult_cell = cell(1,2);
                    
                    % dimension
                    xk = symvar(ch(ii));
                    dim = 1; %for constant term
                    for k=1:length(x)
                        if isequaln(xk,x(k))
                            dim = k;
                            break
                        end
                    end
                    mult_cell{1,2} = dim;
                    
                    % function
                    c_func = matlabFunction(ch(ii),'Vars',{x(dim)});
                    mult_cell{1,1} = c_func;
                    
                    add_cell{ii} = mult_cell;
                    
                else
                    % case: several variable terms multiplied
                    
                    cch = children(ch(ii)); %mult. terms
                    
                    %create cell for mult. terms
                    mult_cell = cell(length(cch),2);
                    for jj = 1:length(cch)
                        
                        % dimension
                        xk = symvar(cch(jj));
                        for k=1:length(x)
                            if isequaln(xk,x(k))
                                dim = k;
                                break
                            end
                        end
                        mult_cell{jj,2} = dim;
                        
                        % function
                        c_func = matlabFunction(cch(jj),'Vars',{x(dim)});
                        mult_cell{jj,1} = c_func;
                        
                    end
                    
                    add_cell{ii} = mult_cell;
                end
                
            end
            
            fcell{i,j} = add_cell;
            
        else
            % case: just one addition term with zero to several variables.
            
            if length(symvar(c_sym)) < 2
                % one variable or constant term
                
                mult_cell = cell(1,2);
                
                % dimension
                xk = symvar(c_sym);
                dim = 1; %for constant term
                for k=1:length(x)
                    if isequaln(xk,x(k))
                        dim = k;
                        break
                    end
                end
                mult_cell{1,2} = dim;
                
                % function
                c_func = matlabFunction(c_sym,'Vars',{x(dim)});
                mult_cell{1,1} = c_func;
                
            else
                % several variables
                
                % create cell corresponding to mult. terms
                mult_cell = cell(length(ch),2);
                for jj = 1:length(ch)
                    
                    % dimension
                    xk = symvar(ch(jj));
                    for k=1:length(x)
                        if isequaln(xk,x(k))
                            dim = k;
                            break
                        end
                    end
                    mult_cell{jj,2} = dim;
                    
                    % function
                    c_func = matlabFunction(ch(jj),'Vars',{x(dim)});
                    mult_cell{jj,1} = c_func;
                    
                end
                
            end
            
            add_cell = cell(1,1);
            add_cell{1} = mult_cell;
            
            fcell{i,j} = add_cell;
            
        end
    end
    
end

end




