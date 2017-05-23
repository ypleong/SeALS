function [op] = bc_neumann(op, dim, bsca, n, fd1dim)
% BC_DIRICHLET Set up Neumann boundary conditions boundary conditions on
% operator. More precisely:
%   Modifies an operator so that it doesn't affect the boundaries of its
%   domain. The operator acts as "pure" differentiation on the boundaries
%   (in the relevant dimension respectively) times the boundary scaling.
%   Common boundary points will have the scaling of the lowest dimension,
%   see also MAKEBC.
% Input:
%   op - the operator as ktensor.
%   dim - the dimension along which to modify the operator.
%   n - the number of grid points in each dimension.
%   fd1dim - the differentiation matrix in the corresponding dimension.
% Output:
%   op - the adjusted operator.
%
% OBS: Neumann conditions has not been tested in any example yet.
%
% See also MAKEBCOP.

% Elis Stefansson, Aug 2015

d = length(n);

%% Construct identity on the boundary without interfering lower dimensions
U = cell(1,d);

% don't interfer with boundary for lower dimensions
for i = 1:(dim-1)
    ni = n(i);
    idcut = eye(ni);
    idcut(1,:) = zeros(1,ni);
    idcut(end,:) = zeros(1,ni);
    U{i} = idcut(:);
end

% the boundary for dim
ndim = n(dim);
dim_matrix = zeros(ndim,ndim);

% scaling incorporating the magnitude from the finite difference scheme fd1
fd1_max = max(abs(fd1dim(1,:)));
dim_matrix(1,:) = (bsca(1)/fd1_max)*fd1dim(1,:);
fd1_max = max(abs(fd1dim(ndim,:)));
dim_matrix(ndim,:) = (bsca(2)/fd1_max)*fd1dim(ndim,:);

U{dim} = dim_matrix(:);

for i = (dim+1):d
    ni = n(i);
    id = eye(ni);
    U{i} = id(:);
end

ident = ktensor(U);

%% Eliminate whatever is happening on the boundary currently
sr = ncomponents(op);
dimU = op.U{dim};
for i=1:sr
    el = reshape(dimU(:,i),ndim,ndim);
    el(1,:) = zeros(1,ndim);
    el(ndim,:) = zeros(1,ndim);
    dimU(:,i) = el(:);
end
op.U{dim} = dimU;

%% Combine above
op = op + ident;
op = fixsigns(arrange(op));
end

