function [B, s] = three_way_array_maker(A)
% This function takes a STT representation of a tensor A and saves it to a
% 3 way array and a list of numbers for the s values

% Save lambdas
s = A.lambda;

% Get the rank of the matrix, etc.
D = ndims(A);
N2 = size(A,1); % Note, since saving to an array, all dir have the same # of sample
rank = length(s);

% Initialize the output array
B = zeros(rank, D, N2);

% Fill up the array
for d = 1:D
    B(:,d,:) = A.U{d}';
end


