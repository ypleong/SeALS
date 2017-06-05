function [Y,err,svCell] = tenid(X,tol,k0,kmax,nrmtype,delta,Xnorm,vb)
% Tensor ID.
%
% Input:
%  X = input tensor
%  tol = stop when  ||X-Y|| < tol*||X||
%  k0 = start with projections of size 2^k0
%  kmax = max projection size is 2^kmax
%  nrmtype = 'frob' or 'snorm'
%  Xnorm = norm of X (should match nrmtype)
%  vb = verbose flag
%
% Output:
%  Y = tensor of reduced separation rank
%  err = err per number of random projections performed

  if vb
    fprintf('TENSOR ID:\n')
  end

  D = ndims(X);
  N2 = size(X);

  tol_id = 1e-8;   % internal tolerance for skeletonization

  kk = 2.^(k0:kmax);   % vector of number of random projections
  dk = kk - [0,kk(1:(end-1))];   % number of new projections per iteration
  
  % if number of terms is already less than minimum number of projections,
  % return
  % Also turned off by MR
  %if length(X.lambda) <= min(kk)
  %  Y = X;
  %  fprintf('Not performing tensor ID, rank is low enough\n')
  %  return
  %end
  
  err = zeros(length(kk),1);
  Ymat = zeros(kk(end),length(X.lambda));
  icurr = 1;
  
  svCell = cell(length(kk),1);
  
  for i = 1:length(kk)
    
    % if you can't do better, return
    % MR turned this off...
    %if (i>1) && (length(X.lambda) <= kk(i));
    %  fprintf('No improvement from tensor ID, exiting\n')
    %  Y = X;
    %  return
    %end
    
    % Added by MR
    % Print which iteration we wre on
    if vb; fprintf('Computing tenid for kk = %d:\n', kk(i)); end;
    
    % generate gram matrix with random vectors
    if vb; fprintf('  Forming Gram matrix...');  end;
    Ri = normalize(krandn(D,N2,dk(i),1));
    if i==1
      R = Ri;
    else
      R = R + Ri;
    end
    
    Ymat(icurr:(icurr+dk(i)-1),:) = gram2(Ri,X); 
    icurr = icurr+dk(i);
    if vb; fprintf('finished.\n'); end;
    
    % skeletonize the gram matrix and get the corresponding terms
    if vb; fprintf('  Skeletonizing...'); end;
    %[P,ix] = kskel(Ymat(1:(icurr-1),:),norm(Ymat(1:(icurr-1),:))*tol_id);
    [P,ix] = kskel(Ymat(1:(icurr-1),:),tol_id); % note: no need to include norm of Y, I changed kskel
    Y = getterms(X,ix);
    Y = arrange(Y);
    Y.lambda = ones(size(Y.lambda));
    if vb
      fprintf('finished.\n'); 
      fprintf('  ||Ymat-Ymat(:,ix)*P||_2 / ||Ymat||_2 = %e\n', norm(Ymat-Ymat(:,ix)*P)/norm(Ymat));
    end
    
    % form gram matrices for computing coefficients
    if vb;  fprintf('  Constructing linear system...');  end
    Gy = full(gram(Y));
    Gx = full(gram2(Y,X));
    if vb;  fprintf('finished.\n');  end;
    
    % use SVD for pseudo inverse
    if vb; fprintf('  Solving linear system...');  end;
    [U,S,V] = svd(Gy);
    
    % Dave's pseudoinverse. I prefer to truncate the SVD. I'm not sure it
    % matters which way you go about it.
    %rnk = min(size(S));
    %Gycond = S(1,1)/S(rnk,rnk);
    %S = S(1:rnk,1:rnk) + S(1,1) * (100*eps) * eye(rnk);   % should this just be, e.g., 100*eps not sqrt(eps)?
    %Gypinv = V(:,1:rnk) * diag(1./diag(S(1:rnk,1:rnk))) * U(:,1:rnk)';
    %P = Gypinv * Gx;
    
    % Pseudo inverse by MR.
    rnk = length(S(diag(S)/S(1,1) > tol_id));
    Gycond = S(1,1)/S(rnk,rnk);
    S = S(1:rnk,1:rnk);
    Gypinv = V(:,1:rnk) * diag(1./diag(S(1:rnk,1:rnk))) * U(:,1:rnk)';
    P = Gypinv * Gx;
    
    % Print the condition numbers
    if vb; fprintf('finished.\n'); end
    if vb; fprintf('    (cond(Gy) = %e)\n', Gycond); end;
    if vb; fprintf('    (cond(Gx) = %e)\n', cond(Gx)); end;
    
    % update the coefficients
    if vb; fprintf('  Computing new coefficients...'); end;
    Y.lambda = Y.lambda .* sum(P,2);
    Y = poswts(Y);
    if vb; fprintf(' finished.\n'); end;

    % Condition the representation to avoid disparate svalues
    % Dave: since I use als after I apply the tenid, I don't require this
    % step and comment it out.
    %if vb; fprintf('  Conditioning via ALS sweep...\n'); end;
    %Y = alsi(Y,length(Y.lambda),0,-Inf,delta,[],'maxit',3,'verbose',vb,'B',Y,'density',1);
    %if vb; fprintf('Finished ALS sweep.\n'); end;
    
    if vb; fprintf('  Computing error...'); end;
    if isequal(nrmtype,'frob')
      err(i) = sqrt(abs(Xnorm^2+fnorm(Y)^2-2*iprod(X,Y)))/Xnorm;
    end
    if isequal(nrmtype,'snorm')
      err(i) = snorma(poswts(Y-X))/Xnorm;
    end
    if vb; fprintf('finished.\n'); end;
    
    if vb; fprintf('kk = %d, NTERMS = %d, TOL = %e, ERR = %e\n\n', kk(i), length(ix), tol, err(i)); end;
    
    % I don't think this is sufficient. Isn't there a requirement on dim?
    if (err(i) < sqrt(2)*tol) && (length(Y.lambda) + 10 <= kk(i)) %added by MR
      err = err(1:i);
      fprintf('Tensor ID converged with %d terms, err = %e\n', length(Y.lambda), err(i));
      return
    end
    
  end
  
end