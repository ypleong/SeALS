clear all;

fprintf('which -all gram\n')
which -all gram
fprintf('which -all gram2\n')
which -all gram2
fprintf('which -all gram2log\n')
which -all gram2log

ok = 1;

d = 10;
n = 1024;
r = 10;
A = krandn(d,n,r,1);
B = krandn(d,n,r,1);

g1 = gram(A);
g2 = gram2(A,A);
g3 = gram2log(A,A);

err1 = norm(g1-g2,'fro')/norm(g1,'fro');
err2 = norm(g1-g3,'fro')/norm(g1,'fro');
fprintf('err1 = %e\n', err1)
fprintf('err2 = %e\n', err2)

g4 = gram(A-B);
g5 = gram2(A-B,A-B);
g6 = gram2log(A-B,A-B);
err3 = norm(g4-g5,'fro')/norm(g4,'fro');
err4 = norm(g4-g6,'fro')/norm(g4,'fro');
fprintf('err3 = %e\n', err3)
fprintf('err4 = %e\n', err4)

if max([err1;err2;err3;err4]) > 1e-14
  ok = 0;  fprintf('FAILED\n');  return;
else
  fprintf('PASSED\n')
end