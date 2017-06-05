clear all;

fprintf('\n***test_getterms***\n')

fprintf('which -all getterms\n')
which -all getterms

ok = 1;

d = 10;
n = 1024;
r = 10;
A = krandn(d,n,r,1);

B = A;
iz = [2;5;6;9;10];
B.lambda(iz) = 0;

inz = 1:r;  inz(iz) = [];
A = getterms(A,inz);

err = fnorm(A-B)/fnorm(B);
fprintf('err = %e\n' ,err)

if err > 1e-7
  ok = 0;  fprintf('FAILED\n');  return;
else
  fprintf('PASSED\n')
end