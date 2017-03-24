function [n,grid,h,region,ord_of_acc,noise_cov,R] = makebasic(bdim,n,region,ord_of_acc,noise_cov,R,G,B)
% MAKEBASIC initalize basic variables. For example, if ord_of_acc is
% specified as a scalar (same order of accuracy in all dimensions) then
% MAKEBASIC converts ord_of_acc into the general form. 
% Inputs:
%   bdim - the dimensions of the hyperrectangle domain.
%   n - the number of grid points in each dimension.
%   region - the dimension of the goal region.
%   ord_of_acc - the order of accuracy for the derivatives.
%   noise_cov - the variance of the noise.
%   R - the control cost matrix (may be scalar).
%   G - the function G in symbolic representation.
%   B - the function B in symbolic representation.
%   See MAIN_PROGRAM for a detailed description.
% Outputs:
%   n - number of grid points in each dimension.
%   grid - the discretization grid.
%   h - the grid step in each dimension.
%   region - the dimensions of the goal region.
%   ord_of_acc - order of accuracy for the derivates.
%   noise_cov - the noise covariance as a matrix.
%   R - the control cost matrix as a matrix.
%
% See also, MAIN_RUN.

% Elis Stefansson, Aug 2015

%% grid

size_bdim = size(bdim);
d = size_bdim(1);

if length(n(:)) == 1
    n = repmat(n,d,1);
end

grid = cell(1,d);
h = ones(1,d);
for i=1:d
    grid{i} = linspace(bdim(i,1),bdim(i,2),n(i))';
    gridi = grid{i};
    h(i) = gridi(2)-gridi(1);
end

%% region

if isempty(region) == 0
    
    % same boundary for all dimensions
    if length(region(:)) == 2
        %check if true dim
        if size(region) ~= [1 2]
            error('wrong size on region')
        end
        region = repmat(region,d,1);
    end
    
    % point case for region, change region to closest point.
    if region(:,1) == region(:,2)
        for i=1:d
            gridi = grid{i};
            [c,index] = min(abs( gridi-region(i,1) ));
            region(i,:) = [gridi(index) gridi(index)];
        end
    end
    
end

%% order of accuracy for derivatives

% same accuracy for all dimensions
if length(ord_of_acc(:)) == 1
    k = ord_of_acc;
    ord_of_acc = [k k; k k];
elseif length(ord_of_acc(:)) == 2
    ord_of_acc = [ord_of_acc; ord_of_acc];
end

%% noise covariance
sizeB = size(B);
if length(noise_cov) == 1
    noise_cov = noise_cov*eye(sizeB(2));
else %then vector for the diagonal elements. 
    noise_cov = diag(noise_cov);
end

%% matrix R
sizeG = size(G);
if length(R) == 1
    R = R*eye(sizeG(2));
end
end