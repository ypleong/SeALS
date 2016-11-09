function res = SRMultM(A,B)
% SRMULTM multiplies two ktensor operators A and B. The result is a new
% ktensor AB.
% Input:
%   A - ktensor operator A.
%   B - ktensor operator B.
% Output:
%   AB - the product of A and B.
%
% Note:
% The matrices of A and B are stored as vectors. For
% example, a matrix of size mxm is stored as a vector with m^2 entries.
%
% SRMULTM is for example used in HESSIAN.

rA = ncomponents(A);
rB = ncomponents(B);
nd = ndims(A);

LL = blockTransposeV2H(A.lambda*B.lambda',1)';
AB = cell(nd,1);
N = round(sqrt(size(B,1)));
for nn = 1:nd
    AB{nn} = reshape(blockTransposeV2H(blockTransposeH2V(reshape(A.U{nn},N,N*rA),N)*reshape(B.U{nn},N,N*rB),N),N*N,[]);
end
res = ktensor(LL,AB);