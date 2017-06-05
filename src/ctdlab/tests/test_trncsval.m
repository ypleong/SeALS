clear all

fprintf('\n***test_trncsval***\n')

fprintf('which -all trncsval\n')
which -all trncsval

ok = 1;

d = 10;
n = 1024;
r = 10;
A = krandn(d,n,r,1);

A = arrange(A);
A.lambda = A.lambda .* exp(-50*linspace(0,1,r).^2)';
A.lambda = A.lambda / A.lambda(1);

B = trncsval(A,1e-10);
err = fnorm(A-B)/fnorm(A);
fprintf('err = %e\n', err)

if err < 1e-7
  fprintf('PASSED\n')
else
  ok = 0;
  fprintf('FAILED\n')
end