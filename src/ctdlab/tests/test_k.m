clear all;

fprintf('\n***test_k***\n')

fprintf('which -all kzeros\n')
which -all kzeros
fprintf('which -all kones\n')
which -all kones
fprintf('which -all knnz\n')
which -all knnz
fprintf('which -all knumel\n')
which -all knumel

ok = 1;

d = 3;
n = 5;
r = 10;

A = full(kzeros(d,n));
Af = zeros([n,n,n]);
res1 = A-Af;
err1 = sqrt(sum(res1(:).^2));
fprintf('err1 = %e\n', err1)

B = full(kones(d,n));
Bf = ones([n,n,n]);
res2 = B-Bf;
err2 = sqrt(sum(res2(:).^2)) / sqrt(sum(Bf(:).^2));
fprintf('err2 = %e\n', err2)

C = kones(d,n);
nnz1 = knnz(C);
nnz2 = knumel(C);
err3 = abs(nnz1-nnz2);
fprintf('err3 = %e\n', err3)

if max(err1,err2) > 1e-15   ||   err3 ~= 0
  ok = 0;  fprintf('FAILED\n');  return;
else
  fprintf('PASSED\n')
end