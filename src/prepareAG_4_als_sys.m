function [AtA, AtG] = prepareAG_4_als_sys(A, G)

G = arrange(G);
A = arrange(A);

nd = ndims(G);
rA = ncomponents(A);
nn = size(G);

AtA = cell(nd,1);
AtG = cell(nd,1);

for d = 1:nd 
    nf = nn(d);
    A_Ud = reshape(A.U{d}*diag((A.lambda).^(1/nd)),nf,nf*rA);
    AtA{d} = A_Ud'*A_Ud;
    AtG{d} = reshape(A_Ud'*G.U{d}*diag((G.lambda).^(1/nd)),nf,[],1);
end
