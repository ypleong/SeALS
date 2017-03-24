function [At] = blockTransposeH2V(A, blksize)

sizeA = size(A);

At = reshape(permute(reshape(A, sizeA(1), blksize, []), [2 1 3]), blksize,[])';  
