function ip=iprod(A,B)
% Inner product of two ktensors in vector format

    C = A.lambda * B.lambda';
    for d = 1:ndims(A)
        C = C .* (A.U{d}'*B.U{d});
    end
    ip = sum(C(:));

end