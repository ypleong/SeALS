clear all

fprintf('\n***test_trncel***\n')

fprintf('which -all trncel\n')
which -all trncel
fprintf('which -all trncel2\n')
which -all trncel2

ok = 1;

d = 10;
n = 1024;
r = 10;
A = krandn(d,n,r,1);
A = (1/fnorm(A)) * A;

tol = 1e-6;

try 
  B = trncel(A,0);
catch
  ok = 0; fprintf('FAILED  :  B = trncel(A,0)\n');  return;
end
try
  B = trncel(A,[]);
catch
  ok = 0; fprintf('FAILED  :  B = trncel(A,[])\n');  return;
end
try
  B = trncel(A,tol);
catch
  ok = 0; fprintf('FAILED  :  B = trncel(A,tol)\n');  return;
end
err = fnorm(A-B)/fnorm(A);
if err > tol
  ok = 0;  fprintf('FAILED  :  err > tol\n');  return;
end

try 
  B = trncel2(A,0);
catch
  ok = 0; fprintf('FAILED  :  B = trncel2(A,0)\n');  return;
end
try
  B = trncel2(A,[]);
catch
  ok = 0; fprintf('FAILED  :  B = trncel2(A,[])\n');  return;
end
try
  B = trncel2(A,tol);
catch
  ok = 0; fprintf('FAILED  :  B = trncel2(A,tol)\n');  return;
end
err = fnorm(A-B)/fnorm(A);
if err > tol
  ok = 0;  fprintf('FAILED  :  err > tol\n');  return;
end

fprintf('PASSED\n')
