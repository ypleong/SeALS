function [F, err, iter, e_list, t_step, illcond, noreduce] = als2(op,varargin)
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

dim = ndims(op);
R = ncomponents(op);
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
err = zeros(maxit,1);
t_step = zeros(maxit,1);
tol = 1000;
noreduce = 0;

e_list = 1;

while (iter < maxit) && (tol > e) && (curr_rank < R)
    
    [P,~,out]= cp_als(op,curr_rank,'init',Pinit,'printitn',0);
    iter = iter + out.iters;
    e_list = [e_list out.iters+1];
    tol = out.err(out.iters);
    illcond = out.ill_cond;
    err(iter-out.iters+1:iter) = out.err;
    t_step(iter-out.iters+1:iter) = out.t_step;
    temp = cell(1,dim);
    for ii = 1:dim
        temp{ii} = matrandnorm(size(P.U{ii},1),1);
        Pinit{ii} = [P.U{ii} temp{ii}];     
    end
    curr_rank = curr_rank + 1;
    
end

t_step = t_step(1:iter);
err = err(1:iter);

if curr_rank == R
    F = op;
    if iter == 0
        err = 0;
    else
        err(iter) = 0;
    end
    disp('cannot reduce rank')
    noreduce = 1;
else
    F = P;
end

if (iter == maxit && err > e)
    disp('Max iterations reached without finding a good solution')
end


