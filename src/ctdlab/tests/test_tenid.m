fprintf('\n***test_tenid***\n')

fprintf('which -all snorma\n')
which -all snorma
fprintf('which -all rofpi1\n')
which -all rofpi1
fprintf('which -all tenid\n')
which -all tenid

ok = 1;

d = 10;
n = 64;
r = 1024;

% generate ctd with exponentially decaying svals
A = krandu(d,n,r,1);
A = normalize(A);
A.lambda = A.lambda / A.lambda(1);
A.lambda = A.lambda .* exp(-1000*linspace(0,1,r).^2)';

% use tensor ID to grab the significant terms, using frobenius error
tol = 1e-7;
k0 = 1;
kmax = 9;
nrmtype = 'frob';
delta = [];
Anorm = fnorm(A);
vb = 0;
fprintf('Target CTD: %d terms above tol\n', length(find(A.lambda>tol)));
fprintf('Running TENID with frobenius norm:\n')
[B,err] = tenid(A,tol,k0,kmax,nrmtype,delta,Anorm,vb);

% now use snorm
fprintf('Running TENID with snorm:\n')
nrmtype = 'snorm';
[B,err] = tenid(A,tol,k0,kmax,nrmtype,delta,Anorm,vb);

