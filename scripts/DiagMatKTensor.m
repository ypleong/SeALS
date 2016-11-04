function [MKTensor_Op] = DiagMatKTensor(MKTensor)
% MATKTENS converts a cell array of ktensor functions to corresponding cell
% array of ktensor operators.
% Input:
%   MKTensor - a cell array with ktensor functions
% Output:
%   MKTensor - corresponding cell array with ktensor operators.
%
% See MAKEOP for an example when it is used.

% Elis Stefansson, Aug 2015

MKTensor_size = size(MKTensor);
MKTensor_Op = cell(MKTensor_size);

for i = 1:MKTensor_size(1)
    for j = 1:MKTensor_size(2)
        MKTensor_Op{i,j} = DiagKTensor(MKTensor{i,j});
    end
end
        