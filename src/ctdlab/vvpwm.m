function c = vvpwm(a,b)
% "Vector-vector pointwise multiplication"

  D = ndims(a);
  ra = length(a.lambda);
  rb = length(b.lambda);

  c.U = cell(1,D);
  c.lambda = zeros(ra*rb,1);

  for d = 1:D
    c.U{d} = zeros(size(a,d),ra*rb);
    ii=1;
    for r = 1:ra
      for s = 1:rb
        c.U{d}(:,ii) = a{d}(:,r).*b{d}(:,s);
        ii=ii+1;
      end
    end
  end
  ii=1;
  for r=1:ra
    for s=1:rb
      c.lambda(ii) = a.lambda(r) * b.lambda(s);
      ii=ii+1;
    end
  end

  c = ktensor(c.lambda,c.U);
  c = arrange(c);

end