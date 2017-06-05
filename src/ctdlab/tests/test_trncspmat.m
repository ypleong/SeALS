clear all

fprintf('\n***test_trncspmat***\n')

fprintf('which -all trncspmat\n')
which -all trncspmat

ok = 1;

n = 512;
delta = 1e-5;

A = sprandn(n,n,.8);

Af = full(A);
Af(abs(Af)<delta) = 0;

B = trncspmat(A,delta);
Bf = full(B);

err = norm(Af-Bf,'fro')/norm(Af,'fro');
fprintf('err = %e\n', err)

if err < 1e-15
  fprintf('PASSED\n')
else
  ok = 0;
  fprintf('FAILED\n')
end