clear all

fprintf('\n***test_poswts***\n')

fprintf('which -all poswts\n')
which -all poswts

ok = 1;

d = 10;
n = 1024;
r = 10;
A = krandn(d,n,r,1);

ix = [1;3;5];
A.lambda(ix) = -A.lambda(ix);

B = poswts(A);

err = fnorm(A-B)/fnorm(A);
fprintf('err = %e\n', err)

if err < 1e-7
  fprintf('PASSED\n')
else
  ok = 0;
  fprintf('FAILED\n')
end