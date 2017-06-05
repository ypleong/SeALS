clear all;

fprintf('\n***test_mvm***\n')

fprintf('which -all keye\n')
which -all keye
fprintf('which -all mvm\n')
which -all mvm

ok = 1;

d = 10;
n = 1024;
r = 10;

a = krandn(d,n,r,1);
II = keye(d,n);

b = mvm(II,a);

err = fnorm(a-b)/fnorm(a);
fprintf('err = %e\n', err)

if err > 1e-7
  ok = 0;  fprintf('FAILED\n');  return;
else
  fprintf('PASSED\n')
end