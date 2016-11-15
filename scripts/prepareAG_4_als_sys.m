function [AtA, AtG] = prepareAG_4_als_sys(A, G)

G = arrange(G);
A = arrange(A);

nd = ndims(G);
rG = ncomponents(G);
rA = ncomponents(A);
nf = size(G,1);

A_U = cell(nd,1);
AtA = zeros(nf*rA,nf*rA,nd);
AtG = zeros(nf,rA*rG,nd);

for d = 1:nd 
    A_U{d} = reshape(A.U{d}*diag((A.lambda).^(1/nd)),nf,nf*rA);
    AtA(:,:,d) = A_U{d}'*A_U{d}; 
    AtG(:,:,d) = reshape(A_U{d}'*G.U{d}*diag((G.lambda).^(1/nd)),nf,[],1); 
end
