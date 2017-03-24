function res = SRMultV(A,F)
% SRMULTV multiplies a ktensor operators A with a ktensor function F. The 
% result is a new ktensor function AF.
% Input:
%   A - ktensor operator.
%   F - ktensor function.
% Output:
%   AB - the product of A and B.
%
% Note:
% The ktensor operator A stores it corresponding matrices as vectors. For
% example, a matrix of size mxm is stored as a vector with m^2 entries.

rA = ncomponents(A);
nd = ndims(F);

LL = blockTransposeV2H(A.lambda*F.lambda',1)';
AU = cell(nd,1);
for nn = 1:nd
    N = size(F,nn);
    AU{nn} = blockTransposeV2H(blockTransposeH2V(reshape(A.U{nn},N,N*rA),N)*F.U{nn},N);
end
res = ktensor(LL,AU);