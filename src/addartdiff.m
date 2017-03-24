function [op] = addartdiff(op,conv,diff,artdiff_option,h,D2,fTens)
% ADDARTDIFF adds artificial diffusion to the operator op.
% Input:
%   op - the operator as ktensor.
%   conv - vector corresponding to which dimensions that have a convection
%   term. conv(i) = 1 if dimension i has convection term, otherwise 0.
%   diff - vector corresponding to which dimensions that have a diffusion
%   term. diff(i) = 1 if dimension i has diffusion term, otherwise 0.
%   artdiff_option - cell with current option for artificial diffusion.
%   h - the grid step in each dimension.
%   D2 - cell with second order differentiation operators. D2{i} is
%   differentiation with respect to dimension i.
%   fTens - f as ktensor cell.
% Outputs:
%   op - the new operator (as ktensor) with added artficial diffusion.
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 2015

%% assign
[artdiff_version,artdiff_dims,artdiff_scale] = deal(artdiff_option{:});

%% version
if artdiff_version == 1
    [opArtdiff] = artdiff_ver1(artdiff_dims,artdiff_scale,conv,diff,D2);
elseif artdiff_version == 2
    [opArtdiff] = artdiff_ver2(artdiff_dims,artdiff_scale,conv,diff,D2,fTens,h);
elseif artdiff_version == 3
    [opArtdiff] = artdiff_ver3(artdiff_dims,artdiff_scale,conv,diff,D2,fTens,h);
else
    error('wrong specifications on the artificial diffusion term')
end

%% add to operator
if isempty(opArtdiff) == 0
    op = op+opArtdiff;
else
    fprintf('OBS: no artificial diffusion term was added')
end

end