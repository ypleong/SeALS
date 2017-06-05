% testGS.m
% This script tests the orthogonality of the outputs of modified Gram-
% Schmidt, doubly modified Gram-schmidt, and Matlab's qr algorithm
% Authors
% Matthew Reynolds

% Create a matrix
sz = 15;
A = randn(sz);
[U,~,V] = svd(A);
S = diag(pi*exp(sqrt(3))*10.^-(0:sz-1));

% norm
nA = S(1,1);

% Set tolerance (relative), for methods that use tolerance
tol = 1e-14;

% print the header
fprintf('\n')
fprintf('Test orthogonality\n')
fprintf('   k(A)        MGS||1-Q^t*Q||   DGS||1-Q^t*Q||  rkGS||1-Q^t*Q||   MAT||1-Q^t*Q||\n')

for i = 1:sz
    % Form the matrix with the condition number 10^(i-1)
    A = U(:,1:i)*S(1:i,1:i)*V(:,1:i)';
    
    % Compute QR factorizations
    [Q,R,piv,k] = dgsqr(A,tol*nA);
    [Qm,Rm,pivm,km] = mgsqr(A,tol*nA);
    [Qrk,Rrk,pivrk,krk] = rank_k_gsqr(A,i); % note: play with your target rank to see how this behaves with an unknown rank.
    [Qmat,Rmat] = qr(A);
    
    % print the results
    fprintf('%4.4e %16.4e %16.4e %16.4e %16.4e\n', S(1,1)/S(i,i), norm(eye(size(Qm,2))-Qm'*Qm),...
        norm(eye(size(Q,2))-Q'*Q), norm(eye(size(Qrk,2))-Qrk'*Qrk), norm(eye(size(Qmat,2))-Qmat'*Qmat))
    
end

% Test the reconstructions
fprintf('\n')
fprintf(' Reconstruction errors, relative \n')
fprintf('   k(A)        MGS||A(:,p)-Q*R||   DGS||A(:,p)-Q*R||  rkGS||A(:,p)-Q*R||   MAT||A-Q*R||\n')
% print the results
    fprintf('%4.4e %16.4e %20.4e %20.4e %16.4e\n', S(1,1)/S(i,i), norm(A(:,pivm)-Qm*Rm)/norm(A(:,pivm)),...
        norm(A(:,piv)-Q*R)/norm(A(:,piv)), norm(A(:,pivrk)-Qrk*Rrk)/norm(A(:,pivrk)), norm(A-Qmat*Rmat)/norm(A))
