function [AtA, AtG] = prepareAG_4_als_sys(A, G)

G = arrange(G);
A = arrange(A);

nd = ndims(G);
rG = ncomponents(G);
rA = ncomponents(A);
nf = size(G,1);

AtA = zeros(nf*rA,nf*rA,nd);
AtG = zeros(nf,rA*rG,nd);

for d = 1:nd 
    A_Ud = reshape(A.U{d}*diag((A.lambda).^(1/nd)),nf,nf*rA);
    AtA(:,:,d) = A_Ud'*A_Ud; 
    AtG(:,:,d) = reshape(A_Ud'*G.U{d}*diag((G.lambda).^(1/nd)),nf,[],1); 
end
