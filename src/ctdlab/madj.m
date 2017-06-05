function A = madj(A)
% Returns the adjoint of an operator.  If the operator is sparse, this is
% quite slow

    N2 = size(A);
    N = sqrt(N2);

    for d = 1:ndims(A)
        for r = 1:length(A.lambda)
            tmp = reshape(A{d}(:,r),[N(d),N(d)]);
            tmp = tmp';
            A{d}(:,r) = reshape(tmp,[N2(d),1]);
        end
    end

end