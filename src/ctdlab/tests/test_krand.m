clear all;

fprintf('\n***test_krand***\n')

fprintf('which -all krandn\n')
which -all krandn
fprintf('which -all krandb\n')
which -all krandb
fprintf('which -all krandu\n')
which -all krandu

ok = 1;

d = 5;
n = 1024;
r = 10;

% not sure what tests to run here ... I guess just make sure they execute

try
  A = krandn(d,n,r,1);
  A2 = krandn(d,n,r,.1);
catch
 ok = 0;  fprintf('FAILED\n');  return;
end
if knnz(A) ~= knumel(A)
  ok = 0;  fprintf('FAILED\n');  return;
end
fprintf('krandn:  p = .1, knnz(A2)/knumel(A2) = %f\n', knnz(A2)/knumel(A2));

try
  A = krandb(d,n,r,1);
  A2 = krandb(d,n,r,.1);
catch
 ok = 0;  fprintf('FAILED\n');  return;
end
if knnz(A) ~= knumel(A)  % A should be dense
  ok = 0;  fprintf('FAILED\n');  return;
end
fprintf('krandb:  p = .1, knnz(A2)/knumel(A2) = %f\n', knnz(A2)/knumel(A2));

try
  A = krandu(d,n,r,1);
  A2 = krandu(d,n,r,.1);
catch
 ok = 0;  fprintf('FAILED\n');  return;
end
if knnz(A) ~= knumel(A)  % A should be dense
  ok = 0;  fprintf('FAILED\n');  return;
end
fprintf('krandu:  p = .1, knnz(A2)/knumel(A2) = %f\n', knnz(A2)/knumel(A2));

fprintf('PASSED\n')