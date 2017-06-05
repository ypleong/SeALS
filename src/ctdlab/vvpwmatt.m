function c = vvpwmatt(a,b)
%function c = vvpwmatt(a,b)
%   This function computes the vector vector pointwise mult more
%   efficiently

% Get the dimensions and rank of each tensor
D = ndims(a);
ra = length(a.lambda);
rb = length(b.lambda);

% initialize the new tensor, with the correct number of dimensions and the
% appropriate 
c.U = cell(1,D);
c.lambda = zeros(ra*rb,1);

% loop through dimension
for d = 1:D
    % initialize the new vector
    c.U{d} = repmat(b{d},[1,ra]);
    
    % form matrix for element-by-element multiplication with c.U{d}
    multMatrix = repmat(a{d}.',[1,rb]).';
    multMatrix = reshape(multMatrix(:),[size(a,d),ra*rb]);
    c.U{d} = c.U{d}.*multMatrix;
end

% parallel version of the loop
% weird, but the shorthand above doesn't work with parfor, i.e. we can't
% use b{d}, only b.U{d}
% tic;
% parfor d = 1:D
% 
%     % initialize the new vector
%     tempCell{d} = repmat(b.U{d},[1,ra]);
%     
%     % form matrix for element-by-element multiplication with c.U{d}
%     multMatrix = repmat(a.U{d}.',[1,rb]).';
%     multMatrix = reshape(multMatrix(:),[size(a,d),ra*rb]);
%     tempCell{d} = tempCell{d}.*multMatrix;
% 
% end
% toc;
%
%c.U = tempCell;

% combine the lambda values
c.lambda = repmat(b.lambda,[ra,1]);
multVect = repmat(a.lambda.',[rb,1]);
multVect = reshape(multVect(:),[ra*rb,1]);
c.lambda = c.lambda.*multVect;

% Form the tensor and output
c = ktensor(c.lambda,c.U);
%c = arrange(c);
end

