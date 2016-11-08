function [] = als_test_code()

load('./test/test_data_VTOL.mat')

err = test_compute_WB(A,G);



end


function [err] = test_compute_WB(A,G)
    
    d = ndims(A);
    rA = ncomponents(A);
    N = sqrt(size(A,1));
    
    WB = cell(rA,rA,d);
    lambda = zeros(rA,rA);
    
    tic
    for kk = 1:d
        A_temp = reshape(A.U{kk},N,N,rA);
        for ii = 1:rA
            for jj = 1:rA 
                lambda(ii,jj) = A.lambda(ii)*A.lambda(jj);
                WB{ii,jj,kk} = lambda(ii,jj)*A_temp(:,:,ii)'*A_temp(:,:,jj);    
            end
        end
    end
    toc
    
    W = als_sys_weights(A,G);
    WB2 = W{1};
    
    err = zeros(rA,rA,d);
    for kk = 1:d
        for ii = 1:rA
            for jj = 1:rA          
                err(ii,jj,kk) = max(max(WB{ii,jj,kk}-WB2((ii-1)*N+(1:N),(jj-1)*N+(1:N),kk)));
            end
        end
    end
end

function [] = test_compute_Wb(b)

end

function [] = test_compute_B(A,F)




end