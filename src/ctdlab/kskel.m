function [P,ix] = kskel(Y,tol)
% Compute the skeletonization of a matrix Y, s.t., Y = Y(:,ix)*P
   
  if issparse(Y)
      Y = full(Y);
  end
  
  % Compute the norm of the input matrix Y
  nY = norm(Y);
  
  % Rank-revealing QR
  %[~,R,ix,k] = mgsqr(Y,tol); % Changed by MR
  [~,R,ix,k] = dgsqr(Y,tol*nY); % Doubly-modified GS to maintain orthogonality
  Pr = eye(size(Y,2));
  Pr = Pr(:,ix);
  ix = ix(1:k);

  [Ru,Rs,Rv] = svd(R(1:k,1:k));
  
  % Dave's pseudo inverse
  %rnk = min(size(Rs));
  %Rs = Rs(1:rnk,1:rnk);
  %Rs = Rs + Rs(1,1) * eps * eye(rnk);
  %Rpinv = Rv(:,1:rnk) * diag(1./diag(Rs)) * Ru(:,1:rnk)';
  %T = Rpinv * R(1:k,(k+1):end);
  
  % My pseudo inverse
  rnk = length(Rs(diag(Rs)/Rs(1,1) > tol));
  Rs = Rs(1:rnk,1:rnk);
  Rpinv = Rv(:,1:rnk) * diag(1./diag(Rs)) * Ru(:,1:rnk)';
  T = Rpinv * R(1:k,(k+1):end);

  P = [eye(k),T];
  P = P * Pr';
   
end