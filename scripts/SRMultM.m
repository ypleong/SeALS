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

%%% Old comments %%%
% A and B are matrices stored with the matrices unrolled into n^2 vectors
% computes AB
%%% Old coments %%%%

rB = ncomponents(B);
rA = ncomponents(A);
nd = ndims(A);

Alambda = A.lambda;
Blambda = B.lambda;

for k = 1:nd
    nA = round(sqrt(length(A.U{k}(:,1))));
    nB = round(sqrt(length(B.U{k}(:,1))));
    
    AU{k} = reshape(A.U{k},nA,nA,rA);
    BU{k} = reshape(B.U{k},nB,nB,rB);
end

first = 1;
for kA = 1:rA
	for kB = 1:rB
		clear T U mat
		U = cell(1,nd);
		for k = 1:nd
			mat = AU{k}(:,:,kA)*BU{k}(:,:,kB);
			U{k} = mat(:);
		end
		T = ktensor(Alambda(kA)*Blambda(kB),U);
		if first
			res = T;
			first = 0;
		else
			res = res + T;
		end
	end
end