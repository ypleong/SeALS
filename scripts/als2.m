function [F, err, iter] = als2(op,varargin)
% ALS code to perform the ALS algorithm as described in Beylkin
% and Mohlenkamp 2005.
% Input:
%   op - tensor to be approximated.
%   optional - optional inputs:
%      e - desired accuracy.
%      maxit - maximum number of iterations.
%      Pinit - inital guess (cell array).
%      if options is empty, ALS will use the default setup maxit = 2000, e
%      = 1e-3 and Pinit = random tensor.
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

debug = 0;

dim = ndims(op);
op_len = length(varargin);

if op_len < 3
    Pinit = cell(1,dim);
    for ii = 1:dim
        Pinit{ii} = matrandnorm(size(op.U{ii},1),1);
    end
else
    Pinit = varargin{3};
end

if op_len < 2
    maxit = 2000;
else
    maxit = varargin{2};
end

if op_len < 1
     e = 1e-4; 
else
    e = varargin{1};
end

curr_rank = 1;
iter = 0;
if debug
    err = zeros(maxit,1);
end
tol = 1000;

while (iter < maxit) && (tol > e)
    
    [P,~,out]= cp_als(op,curr_rank,'init',Pinit,'printitn',0);
    iter = iter + out.iters;
    tol = out.err;
    if debug
        err(iter) = tol;
    end
    temp = cell(1,dim);
    for ii = 1:dim
        temp{ii} = matrandnorm(size(P.U{ii},1),1);
        Pinit{ii} = [P.U{ii} temp{ii}];     
    end
    curr_rank = curr_rank + 1;
    
end

if ~debug
    err = tol;
end

if (iter == maxit && err > e)
    disp('Max iterations reached without finding a good solution')
end

F = P;

