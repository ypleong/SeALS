clear all;

fprintf('\n***test_vvpwm***\n')

fprintf('which -all vvpwm\n')
which -all vvpwm

fprintf('which -all vvpwmatt\n')
which -all vvpwmatt

ok = 1;

d = 10;
n = 1024;
r = 10;

A = kones(d,n);
B = krandu(d,n,r,1);

tic;
C = vvpwm(A,B);
tC = toc;

% Added by Mr
tic;
D = vvpwmatt(A,B);
tD = toc;

err = fnorm(B-C)/fnorm(B);
errD = fnorm(B-D)/fnorm(B);
fprintf('err vvpwm = %e\n', err)
fprintf('err vvpwmatt = %e\n', errD)

fprintf('Time taken by vvpwm = %e\n',tC)
fprintf('Time taken by vvpwmatt = %e\n',tD)


if err > 1e-7 || errD > 1e-7
  ok = 0;  fprintf('FAILED\n');  return;
else
  fprintf('PASSED\n')
end