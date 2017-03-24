function [At] = blockTransposeV2H(A, blksize)

sizeA = size(A);

At = reshape(permute(reshape(A', sizeA(2), blksize, []), [2 1 3]), blksize,[]);  
