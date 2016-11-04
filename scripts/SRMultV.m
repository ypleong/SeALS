function res = SRMultV(A,F)
% SRMULTM multiplies a ktensor operators A with a ktensor function F. The 
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

rF = ncomponents(F);
rA = ncomponents(A);
nd = ndims(F);

% tensor access is slow, so do this at the expense of mem
Alambda = A.lambda;
Flambda = F.lambda;
AU = cell(nd,1);
FU = cell(nd,1);
for k = 1:nd
    n = size(F.U{k},1);
	AU{k} = reshape(A.U{k},n,n,rA);
	FU{k} = F.U{k};
end


first = 1;
for kA = 1:rA
	for kF = 1:rF
		clear T U
		U = cell(1,nd);
		for k = 1:nd
			U{k} = AU{k}(:,:,kA)*FU{k}(:,kF);
		end
		T = ktensor(Alambda(kA)*Flambda(kF),U);
		if first
			res = T;
			first = 0;
		else
			res = res + T;
		end
	end
end