function [tens] = DiagKTensor(vec)
% DIAGKTENSOR Diagonalize a series of vectors in each dimension to a
% diagonal tensor operator
%Input:
%   vec - ktensor
%Output:
%   tens - corresponding ktensor operator

d = ndims(vec);
sr = ncomponents(vec);

% incooporate scale factors into vec, evenly distributed. Then we
% don't need to worry about vec.lambda
for i=1:d
    vec.U{i} = vec.U{i}*diag((vec.lambda.^(1/d)));
end
vec.lambda = ones(sr,1);


Ut = cell(1,d);
for i=1:d
    myDiag = diag(vec.U{i}(:,1));
    Ut{i} = myDiag(:);
end
tens = ktensor(Ut);
 
for k=2:sr
    Ut = cell(1,d);
    for i=1:d
        myDiag = diag(vec.U{i}(:,k));
        Ut{i} = myDiag(:);
    end

    tens = tens + ktensor(Ut);
end

end

