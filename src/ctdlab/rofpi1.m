function [T,err,alpha] = rofpi1(A,T,tol,maxit,vbs)
% Fixed point iteration for rank-1 update with provided initial guess.  
% If any of the last three inputs are empty, default values are assigned.
%
% Input:
%   A = tensor whose rank-1 approximation is sought
%   T = initial guess for iteration
%   tol = stopping tolerance for change in estimated weight
%   maxit = max number of iterations
%   vbs = verbose flag
% Output:
%   T = rank-1 approximation
%   err = change in weight per iteration
%   alpha = value of weight at each iteration
    
  % Defaults
  if isempty(tol)
      tol = 1e-8;
  end
  if isempty(maxit)
      maxit = 20;
  end

  if nargin<5
    vbs = 0;
  end

  % output variables
  err = zeros(1,maxit);
  alpha = zeros(1,maxit);

  for iter = 1:maxit

    if vbs
      fprintf('iter = %d\n', iter)
    end

    for k = 1:ndims(A)

        T{k} = iprodk(A,T,k);
        alpha(iter) = sqrt(sum(T{k}.^2));
        T{k} = T{k} / alpha(iter);
        T.lambda = alpha(iter);

    end

    if iter > 1
        err(iter) = abs(alpha(iter-1)-alpha(iter));
    end
    
    if vbs
      fprintf('iter = %d, alpha = %e, dalpha = %e\n', iter,...
        alpha(iter), err(iter));
    end
    
    if err(iter) < tol
      err = err(1:iter);
      alpha = alpha(1:iter);
      return
    end

  end


end