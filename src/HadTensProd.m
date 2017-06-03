function res = HadTensProd(F,G)
% HadTensProd performs Hadamand tensor products for a ktensor function F 
% and a ktensor function G. The result is a new ktensor function FG.
% Input:
%   F - ktensor function.
%   G - ktensor function.
% Output:
%   FG - the Hadamand product of F and G.

nd = ndims(F);

if nd ~= ndims(G)
    error('The dimensions of the ktensors must be the same')
end

if size(F) ~= size(G)
    error('The sizes of the ktensors must be the same')
end

rF = ncomponents(F);
rG = ncomponents(G);

LL = G.lambda*F.lambda';
FG = cell(nd,1);
for nn = 1:nd
    N = size(F,nn);
    F_U = F.U{nn};
    temp = zeros(N,rF*rG);
    for ii = 1:rF
        temp(:,((ii-1)*rG)+(1:rG)) = bsxfun(@times, G.U{nn}, F_U(:,ii));
    end
    FG{nn} = temp;
end
res = ktensor(LL(:),FG);

res = fixsigns(arrange(res));