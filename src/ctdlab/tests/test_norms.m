clear all;

fprintf('\n***test_norms***\n')

% make sure you're calling the functions you think you are
fprintf('which -all krandn\n')
which -all krandn
fprintf('which -all fnorm\n')
which -all fnorm

ok = 1;

d = 10;
n = 1024;
r = 3;
A = krandn(d,n,r,1);

norm1 = norm(A);
norm2 = fnorm(A);
err = abs(norm1-norm2)/norm2;

fprintf('Comparing "fnorm" with tensor toolbox "norm" for small ctd:\n')
fprintf('rel err = %e\n', abs(norm1-norm2)/norm2);

if err > 1e-14;  
  ok = 0;
  fprintf('FAILED\n')
  return
else
  fprintf('PASSED\n')
end