function C = mmm(A,B,del)
% Multiply two operators in canonical format, performing multiplication for
% only those terms whose svalues are above a threshold.

  if nargin < 3
    del = 0;
  end

  D = ndims(A);
  N2 = size(A);
  N = sqrt(N2);

  % Compute the new svalues before multiplication / normalization
  lambda = A.lambda * B.lambda';
  [ia,ib,~] = find(abs(lambda) > (max(abs(lambda(:)))*del));
  
  % loop over all terms and perform multiplication
  C.lambda = zeros(length(ia),1);
  C.U = cell(1,ndims(A));
  for d = 1:D
    C.U{d} = sparse(N2(d),length(C.lambda));
    ii = 1;
    for r = 1:length(ia)
      Atmp = reshape(A{d}(:,ia(r)),[N(d),N(d)]);
      if d==1
        C.lambda(ii) = lambda(ia(r),ib(r));
      end
      Btmp = reshape(B{d}(:,ib(r)),[N(d),N(d)]);
      Ctmp = Atmp * Btmp;     
      lamtmp = sqrt(sum(Ctmp(:).^2));
      C.lambda(ii) = C.lambda(ii) * lamtmp;
      C.U{d}(:,ii) = reshape(Ctmp,[N2(d),1]) / lamtmp;
      ii = ii+1;
    end
  end

  C = ktensor(C.lambda,C.U);
    
end