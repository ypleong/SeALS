function nrm = fnorm(A)
% Frobenius norm of ktensor.
%
% Required input:
%   A = ktensor
%
% Output:
%   nrm = Frobenius norm

  nrm = iprod(A,A);
  nrm = full(sqrt(abs(nrm)));
    
end