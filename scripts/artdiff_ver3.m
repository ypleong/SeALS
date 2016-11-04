function [opArtdiff] = artdiff_ver3(artdiff_dims,artdiff_scale,conv,diff,D2,fTens,h)
% ARTDIFF_VER1 is the third version for adding artificial diffusion.
% This version does scale with the magnitude of the components of f 
% (seen as the the convection). See the paper for details.
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

for i=1:d
    if norm(fTens{i}) ~= 0
        
        f_sqr_dim_i = SRMultV(DiagKTensor(fTens{i}),fTens{i});
        f_sqr_dim_i = arrange(f_sqr_dim_i);
        
        % incooporate scale factor in f_sqr_dim_i, evenly distributed
        for k = 1:d
            f_sqr_dim_i.U{k} = f_sqr_dim_i.U{k}*diag((f_sqr_dim_i.lambda).^(1/d));
        end
        f_sqr_dim_i.lambda = ones(ncomponents(f_sqr_dim_i),1);
        
        % approximate |f_i| by sqrt of dominant term in |f_i|^2
        U = cell(1,d);
        for k = 1:d
            U{k} = sqrt( f_sqr_dim_i.U{k}(:,1) );
        end
        f_abs_dim_i = ktensor(U);
        
        % add to artificial diffusion term
        if strcmp(artdiff_dims,'all')
            
            if isempty(opArtdiff) == 1
                opArtdiff = SRMultM(DiagKTensor(f_abs_dim_i),D2{i});
            else
                opArtdiff = opArtdiff+SRMultM(DiagKTensor(f_abs_dim_i),D2{i});
            end
            
        elseif strcmp(artdiff_dims,'needed')
            
            if conv(i) == 1 && diff(i) == 0
                if isempty(opArtdiff) == 1
                    opArtdiff = SRMultM(DiagKTensor(f_abs_dim_i),D2{i});
                else
                    opArtdiff = opArtdiff+SRMultM(DiagKTensor(f_abs_dim_i),D2{i});
                end
            end
            
        else
            error('wrong specifications on the artificial diffusion term');
        end
        
    end
end

if isempty(opArtdiff) == 0
    opArtdiff = (0.5*artdiff_scale*h_max)*opArtdiff;
end