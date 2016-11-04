function [opArtdiff] = artdiff_ver1(artdiff_dims,artdiff_scale,conv,diff,D2)
% ARTDIFF_VER1 is the first version for adding artificial diffusion.
% This version does not scale with the magnitude of f (seen as the 
% convection). See the paper for details.
% Input:
%   artdiff_dims - 'all' if all dimension shall have artificial diffusion
%   and 'needed' if just dimensions that has a convection term and no 
%   corresponding diffusion term shall have artificial diffusion.
%   artdiff_scale - the artificial diffusion scaling parameter, alpha in
%   the paper.
%   conv - vector corresponding to which dimensions that have a convection
%   term. conv(i) = 1 if dimension i has convection term, otherwise 0.
%   diff - vector corresponding to which dimensions that have a diffusion
%   term. diff(i) = 1 if dimension i has diffusion term, otherwise 0.
%   D2 - cell with second order differentiation operators. D2{i} is
%   differentiation with respect to dimension i.
% Output:
%   opArtdiff - the artificial diffusion term as ktensor.
%
% See also ADDARTDIFF.

% Elis Stefansson, Aug 2015

d = length(conv);

if strcmp(artdiff_dims,'all')
    
    opArtdiff =  D2{1};
    for i=2:d
        opArtdiff = opArtdiff+D2{i};
    end
    
    opArtdiff = artdiff_scale*opArtdiff;
    
elseif strcmp(artdiff_dims,'needed')
    
    opArtdiff = [];
    for i=1:d
        if conv(i) == 1 && diff(i) == 0
            if isempty(opArtdiff) == 1
                opArtdiff = D2{i};
            else
                opArtdiff = opArtdiff+D2{i};
            end
        end
    end
    
    if isempty(opArtdiff) == 0
        opArtdiff = artdiff_scale*opArtdiff;
    end
    
else
    error('wrong specifications on the artificial diffusion term')
end







