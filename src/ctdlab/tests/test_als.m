fprintf('\n***test_als***\n')

fprintf('which -all als\n')
which -all als
fprintf('which -all alsi\n')
which -all alsi

N=16;

x1=2*pi*linspace(0,1,N)';
x2=x1;
x3=x1;
x4=x1;
s1=sin(x1);
s2=sin(x2);
s3=sin(x3);
s4=sin(x4);
c1=cos(x1);
c2=cos(x2);
c3=cos(x3);
c4=cos(x4);

s=ones(8,1);

U=cell(1,4);
U{1}=[s1, -s1, c1, -c1, c1, c1, -s1, -s1];
U{2}=[c2, c2, s2, s2, c2, c2, s2, s2];
U{3}=[c3, s3, c3, s3, s3, c3, s3, c3];
U{4}=[c4, s4, c4, s4, c4, s4, c4, s4];

A = ktensor(s,U);

% put als code here
tol = 1e-6;
stucktol = 1e-8;
delta = [];
r0 = 1;
rmax = []; 
verbose = 0;
maxit = 250;
errtype = 'rel';
Anorm = fnorm(A);
A_guess = [];
[B,err,Anorm,ext] = als(A,tol,stucktol,delta,r0,rmax,'verbose',verbose,...
  'maxit',maxit,'errtype',errtype,'Anorm',Anorm,'B',A_guess);
fprintf('ALS: target tensor has rank = %d, optimal reduced srank = %d\n', length(A.lambda), ndims(A))
fprintf('     approximation err = %e, srank = %d\n', err(end), length(B.lambda));