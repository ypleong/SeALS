function [ X ] = mattNormalize( X )
%function [ X ] = mattNormalize( X )
%   This function normalizes using a different formula: sqrt(N)/sqrt(sum(u^2))

%% Normalize matrices
for r = 1:length(X.lambda)
    for n = 1:ndims(X)
        tmp = norm(X.u{n}(:,r),2)/sqrt(length(X.u{n}(:,r)));
        if (tmp > 0)            
            X.u{n}(:,r) = X.u{n}(:,r) / tmp;
        end
        X.lambda(r) = X.lambda(r) * tmp;        
    end
end

%% Ensure that the lambda values are positive
idx = find(X.lambda < 0);
X.u{1}(:,idx) = -1 * X.u{1}(:,idx);
X.lambda(idx) = -1 * X.lambda(idx);



end

