function [Q,R,piv,k] = mgsqr(A,tol)
% Pivoted QR factorization via modified Gram-Schmidt.

  [m,n] = size(A);
  maxrnk = min(m,n);
  Q = zeros(m,maxrnk);
  R = zeros(maxrnk,n);
  
  piv = 1:n;
  c = sum(A.^2,1);
  [cmax,imax] = max(c);
  cmax = cmax(1);  imax = imax(1); % In case two equal col norms

  for k = 1:maxrnk
    
    % perform pivot
    tmp = piv(k);  piv(k) = piv(imax);  piv(imax) = tmp;
    tmp = A(:,k);  A(:,k) = A(:,imax);  A(:,imax) = tmp;
    tmp = R(:,k);  R(:,k) = R(:,imax);  R(:,imax) = tmp;
    tmp = c(k);    c(k) = c(imax);      c(imax) = tmp;
    
    % gram-schmidt sweep
    R(k,k) = norm(A(:,k));
    Q(:,k) = A(:,k) / R(k,k);
    R(k,(k+1):n) = Q(:,k)' * A(:,(k+1):n);
    
    for j = (k+1):n
      A(:,j) = A(:,j) - R(k,j) * Q(:,k);
    end
    
    % find next pivot 
    if k < maxrnk
      %c((k+1):n) = c((k+1):n) - R(k,(k+1):n).^2; % I had problems getting
      %this to work right. As it is written now, you lose 1/2 of the digits
      
      % Slower, but it works
      c((k+1):n) = sqrt(sum(A((k+1):m,(k+1):n).^2,1));
      
      % Find the max column norm
      [cmax,imax] = max(c((k+1):n)); 
      cmax = cmax(1);  imax = imax(1);  imax = k + imax(1);
      
    end
    
    % check for convergence
    if abs(cmax) < tol
      Q = Q(:,1:k);
      R = R(1:k,:);
      return
    end
    
    
  end

end