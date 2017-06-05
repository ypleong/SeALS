clear all;

fprintf('\n***test_mmm***\n')

fprintf('which -all mmm\n')
which -all mmm
fprintf('which -all madj\n')
which -all madj

ok = 1;

d = 10;
n = 32;
n2 = n^2;
r = 10;

% test multiplication from both sides by identity
A = krandn(d,n2,r,1);
II = keye(d,n);

B1 = mmm(II,A);
B2 = mmm(A,II);

err1 = fnorm(A-B1)/fnorm(A);
fprintf('err1 = %e\n', err1)
err2 = fnorm(A-B2)/fnorm(A);
fprintf('err2 = %e\n', err2)
err3 = fnorm(B1-B2)/fnorm(B1);
fprintf('err3 = %e\n', err3);

% test multiplication by both sides of a symmetric operator
Asym = mmm(madj(A),A);
Asym2 = madj(Asym);
err4 = fnorm(Asym-Asym2)/fnorm(Asym);
fprintf('err4 = %e\n', err4);

% test multiplication with truncation
A = krandn(d,n2,r,1);
A = normalize(A);
A.lambda = A.lambda .* exp(-50*linspace(0,1,r).^2)';
B = mmm(A,A);
B2 = mmm(A,A,1e-8);
err5 = fnorm(B-B2) / fnorm(B);
fprintf('err5 = %e, rank(B) = %d, rank(B2) = %d\n', err5,length(B.lambda),length(B2.lambda))

if max([err1,err2,err3,err4,err5]) > 1e-7
  ok = 0;  fprintf('FAILED\n');  return;
else
  fprintf('PASSED\n')
end
