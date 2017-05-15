function [x, DM] = chebdifn(N, M, xmin, xmax, imin, imax)

% Without (a) and (b)

%  The function [x, DM] =  chebdif(N,M) computes the differentiation 
%  matrices D1, D2, ..., DM on Chebyshev nodes. 
% 
%  Input:
%  N:        Size of differentiation matrix.     
%  M:        Number of derivatives required (integer).
%  xmin:     lower bound of the domain
%  xmax:     upper bound of the domain
%  imin:     lower bound of the interior domain
%  imax:     upper bound of the interior domain
%  Note:     0 < M <= N-1.
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
%  The code implements two strategies for enhanced 
%  accuracy suggested by W. Don and S. Solomonoff in 
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric 
%  identities to avoid the computation of differences 
%  x(k)-x(j) and (b) the use of the "flipping trick"
%  which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
%  Note added May 2003:  It may, in fact, be slightly better not to
%  implement the strategies (a) and (b).   Please consult the following
%  paper for details:   "Spectral Differencing with a Twist", by
%  R. Baltensperger and M.R. Trummer, to appear in SIAM J. Sci. Comp. 

%  J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by 
%  JACW, May 2003.

if nargin < 3
    xmin = -1;
    xmax = 1;
    imin = [];
    imax = [];
elseif nargin < 5
    imin = [];
    imax = [];
elseif nargin == 6
    if N < 5
        error('N must be greater than 4')
    end
end
    
     if isempty(imin) && isempty(imax)
         k = (0:N-1)';                        % Compute theta vector.
        th = k*pi/(N-1);
        
        x = cos(th); % Compute Chebyshev points.
        x = ((xmax-xmin)*x+xmax+xmin)/2;
        
        DM = computeD(x, (0:N-1)', N, M);
       
     else
         imid = (imax+imin)/2;
         Nb = floor((imid-xmin)*N/(xmax-xmin));
         
         k1 = (0:Nb)';
         k2 = (0:N-Nb-1)';
         th1 = k1*pi/(Nb);
         th2 = k2*pi/(N-Nb-1);
         
         x1 = cos(th1); % Compute Chebyshev points.
         x2 = cos(th2); % Compute Chebyshev points.
         
         dx = abs(x1(end-1) - x1(end))*(xmax-imid)/2;
         imin = imid-dx;
         imax = imid+dx;
         x1(1) = [];
         x2(end) = [];
               
         x = [((xmax-imax)*x2+imax+xmax)/2; imid; ((imin-xmin)*x1+imin+xmin)/2;];
         
         DM = computeD(x, (0:N-1)', N, M);
         
     end
     
     
end

function [DM] = computeD(x, k, N, M)

     I = eye(N);                          % Identity matrix.     
     L = logical(I);                      % Logical identity matrix.

     T = repmat(x,1,N);                
    DX = T-T';           
 DX(L) = ones(N,1);                       % Put 1's on the main diagonal of DX.

     C = toeplitz((-1).^k);               % C is the matrix with 
C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2;     % entries c(k)/c(j)
C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;

     Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))  
  Z(L) = zeros(N,1);                      % with zeros on the diagonal.

     D = eye(N);                          % D contains diff. matrices.

DM = zeros(N,N,M);

for ell = 1:M
    D = ell*Z.*(C.*repmat(diag(D),1,N) - D); % Off-diagonals
    D(L) = -sum(D,2);                            % Correct main diagonal of D
    DM(:,:,ell) = D;                                   % Store current D in DM
end
end