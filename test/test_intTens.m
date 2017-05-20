function [] = test_intTens()

addpath('../src/')

U{1} = [4 5 6]';
U{2} = [10 20 30]';
tt = ktensor(U);

grid{1} = 1:3;
grid{2} = 1:3;

ww{1} = ones(3,1);
ww{2} = ones(3,1);

assert(intTens(tt) == 400)
assert(intTens(tt,[1 2]) == 400)
assert(intTens(tt,[1 2],grid) == 400)
assert(intTens(tt,[],grid) == 400)
assert(intTens(tt,[1 2],grid,ww) == 400)
assert(intTens(tt,[],[],ww) == 400)
assert(all(double(intTens(tt,1,grid,ww)) == [100 200 300]))
assert(all(double(intTens(tt,2,grid,ww)) == [160 200 240]'))

ww{1} = [1 2 3]';
ww{2} = [4 5 6]';

assert(intTens(tt,[],[],ww) == 4410)
assert(all(double(intTens(tt,1,[],ww)) == [210 420 630]))
assert(all(double(intTens(tt,2,[],ww)) == [840 1050 1260]'))

U{1} = [4 5 6; 4 5 6]';
U{2} = [10 20 30; 10 20 30]';
tt = ktensor(U);

assert(intTens(tt,[],[],ww) == trapz(trapz(double(tt).*repmat(ww{1},1,3).*repmat(ww{2}',3,1))))
assert(intTens(tt,[1 2],[],ww) == trapz(trapz(double(tt).*repmat(ww{1},1,3).*repmat(ww{2}',3,1))))
assert(all(double(intTens(tt,1,[],ww)) == trapz(double(tt).*repmat(ww{1},1,3))))
assert(all(double(intTens(tt,2,[],ww)) == trapz(double(tt).*repmat(ww{2}',3,1),2)))

lambda = [2 6]';
U{1} = [4 5 6; 1 2 3]';
U{2} = [10 20 30; 4 5 6]';
tt = ktensor(lambda,U);

assert(intTens(tt,[],[],ww) == trapz(trapz(double(tt).*repmat(ww{1},1,3).*repmat(ww{2}',3,1))))
assert(intTens(tt,[1 2],[],ww) == trapz(trapz(double(tt).*repmat(ww{1},1,3).*repmat(ww{2}',3,1))))
assert(all(double(intTens(tt,1,[],ww)) == trapz(double(tt).*repmat(ww{1},1,3))))
assert(all(double(intTens(tt,2,[],ww)) == trapz(double(tt).*repmat(ww{2}',3,1),2)))

d = 2;
r = 9;
for ii = 1:d
    n(ii) = randi(100);
    U{ii} = matrandnorm(n(ii),r);
    ww{ii} = randn(n(ii),1);
    grid{ii} = sort([-5; -5+rand(n(ii)-2,1)/10; 5]);
end
tt = ktensor(randn(r,1),U);

fulltt = double(tt);
assert(abs(intTens(tt) - trapz(trapz(fulltt))) <= 100*eps)
assert(abs(intTens(tt,[1 2],[]) - trapz(trapz(fulltt))) <= 100*eps)
assert(all(abs(double(intTens(tt,1)) - trapz(fulltt))<= 100*eps))
assert(all(abs(double(intTens(tt,2)) - trapz(fulltt,2))<= 100*eps))
assert(abs(intTens(tt,[],[],ww) - trapz(trapz(fulltt.*(ww{1}*ww{2}')))) <= 100*eps)
assert(abs(intTens(tt,[1 2],[],ww) - trapz(trapz(fulltt.*(ww{1}*ww{2}')))) <= 100*eps)
assert(all(abs(double(intTens(tt,1,[],ww)) - trapz(fulltt.*repmat(ww{1},1,n(2))))<= 100*eps))
assert(all(abs(double(intTens(tt,2,[],ww)) - trapz(fulltt.*repmat(ww{2}',n(1),1),2))<= 100*eps))
assert(abs(intTens(tt,[],grid,ww) -  trapz(grid{2},trapz(grid{1},fulltt.*(ww{1}*ww{2}')))) <= 100*eps)
assert(abs(intTens(tt,[1 2],grid,ww) - trapz(grid{2},trapz(grid{1},fulltt.*(ww{1}*ww{2}'))))<= 100*eps)
assert(all(abs(double(intTens(tt,1,grid,ww)) - trapz(grid{1},fulltt.*repmat(ww{1},1,n(2))))<= 10*eps))
assert(all(abs(double(intTens(tt,2,grid,ww)) - trapz(grid{2},fulltt.*repmat(ww{2}',n(1),1),2))<= 10*eps))

end