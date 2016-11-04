function [C] = TMatM(A,B)
% TMATM Tensor Matrix Multiply. Multiplies two cell array with ktensor
% operators as elements as if the cell arrays where matrices.
% Inputs:
%   A - a cell array with ktensor operators.
%   B - a cell array with ktensor operators.
% Outputs:
%   C - C=AB, that is corresponding cell array with ktensor operators.
%
% See MAKEOP for an example when it is used.

adim = size(A);
bdim = size(B);

if adim(2) ~= bdim(1)
    error('Cell-matrices must be correctly sized');
end

C = cell(adim(1),bdim(2));

for i=1:adim(1)
    for j=1:bdim(2)
        
        C{i,j} = SRMultM(A{i,1},B{1,j});
        for k=2:adim(2)
            C{i,j} = C{i,j} + SRMultM(A{i,k},B{k,j});
        end
        
    end
end

end

