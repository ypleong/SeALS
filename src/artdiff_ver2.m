function [opArtdiff] =  artdiff_ver2(artdiff_dims,artdiff_scale,conv,diff,D2,fTens,h)
% ARTDIFF_VER1 is the second version for adding artificial diffusion.
% This version does scale with the magnitude of f (seen as the convection).
% See the paper for details.
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
%   h - the grid step in each dimension.
% Output:
%   opArtdiff - the artificial diffusion term as ktensor.
%
% See also ADDARTDIFF.

% Elis Stefansson, Aug 2015

d = length(conv);
h_max = max(h);
opArtdiff = [];

% calculate |f|^2
f_sqr = [];
for i=1:d
    if norm(fTens{i}) ~= 0
        
        if isempty(f_sqr) == 1
            f_sqr = SRMultV(DiagKTensor(fTens{i}),fTens{i});
        else
            f_sqr = f_sqr + SRMultV(DiagKTensor(fTens{i}),fTens{i});
        end
        
    end
end

if isempty(f_sqr) == 0
    
    f_sqr = arrange(f_sqr);
    
    % incooporate scale factor in f_sqr, evenly distributed
    for i = 1:d
        f_sqr.U{i} = f_sqr.U{i}*diag((f_sqr.lambda).^(1/d));
    end
    f_sqr.lambda = ones(ncomponents(f_sqr),1);
    
    % approximate |f| by sqrt of dominant term in |f|^2
    U = cell(1,d);
    for i = 1:d
        U{i} = sqrt( f_sqr.U{i}(:,1) );
    end
    f_abs = ktensor(U);
    
    if strcmp(artdiff_dims,'all')
        
        for i=1:d
            if isempty(opArtdiff) == 1
                opArtdiff = SRMultM(DiagKTensor(f_abs),D2{i});
            else
                opArtdiff = opArtdiff+SRMultM(DiagKTensor(f_abs),D2{i});
            end
        end
        
        if isempty(opArtdiff) == 0
            opArtdiff = (0.5*artdiff_scale*h_max)*opArtdiff;
        end
        
    elseif strcmp(artdiff_dims,'needed')
        
        for i=1:d
            if conv(i) == 1 && diff(i) == 0
                if isempty(opArtdiff) == 1
                    opArtdiff = SRMultM(DiagKTensor(f_abs),D2{i});
                else
                    opArtdiff = opArtdiff+SRMultM(DiagKTensor(f_abs),D2{i});
                end
            end
        end
        
        if isempty(opArtdiff) == 0
            opArtdiff = (0.5*artdiff_scale*h_max)*opArtdiff;
        end
        
    else
        error('wrong specifications on the artificial diffusion term')
    end
    
end

