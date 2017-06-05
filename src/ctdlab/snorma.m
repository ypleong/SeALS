function [nrm,T] = snorma(A,tol,maxit)
  % Estimate the snorm of tensor A using top singular vectors of factor
  % matrices.
  % Inputs:
  %  A = ktensor
  %  tol = stopping tolerance (defaults to 1e-10)
  %  maxit = max iterations (defaults to 50)
  % Outputs:
  %  nrm = s-value of rank-1 approximation (snorm)
  %  T = corresponding rank-1 approximation
  
  if nargin <3
    maxit = 25;
  end
  if nargin < 2
    tol = 1e-8;
  end
  
  T = cell(ndims(A),1);
  for d = 1:ndims(A)
    %[T{d},~,~,~] = powm(A{d},tol,maxit);
    %T{d} = sparse(T{d});
    T{d} = sum(A{d}*diag(A.lambda),2);
    T{d} = T{d} / norm(T{d});
  end
  T = ktensor(1,T);
  
  [T,~,alpha] = rofpi1(A,T,tol,maxit);
  nrm = alpha(end);

end