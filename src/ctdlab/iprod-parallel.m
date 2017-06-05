function ip=iprod(A,B)
% Inner product of two ktensors in vector format
% This version uses a parallel version of a for loop

nd = ndims(A);

parfor d = 1:nd
    mult_cell{d} = A.U{d}'*B.U{d};   
end

C = A.lambda * B.lambda';
%parfor d = 1:ndims(A)
%    C = C .* mult_cell{d};
%end

Cdims = size(C);

for d = 1:ndims(A)
   for i = 1:Cdims(1)
      mc = mult_cell{d}(i,:);
      ci = C(i,:);
      parfor j = 1:Cdims(2) 
        ci(j) = ci(j)*mc(j);
      end
      C(i,:) = ci;
   end
end


% Sum to get the inner product
ip = sum(C(:));

%     C = A.lambda * B.lambda';
%     for d = 1:ndims(A)
%         C = C .* (A.U{d}'*B.U{d});
%     end
%     ip = sum(C(:));

end