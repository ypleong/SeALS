function [] = als_test_code()

% load('./test/test_data_smooth2D.mat')
% load('./test/test_data_VTOL.mat')

d = 4; % ndims(A);
F = cell(1,d);
N = 150; % size(G,1);
rF = 3;
rA = 8;
rG = 4;
for ii = 1:d
    F{ii} = matrandnorm(N,rF);
    A{ii} = matrandnorm(N*N,rA);
    A2{ii} = matrandnorm(N*N,rA+2);
%     G{ii} = matrandnorm(N,rG);
end

At = ktensor(A);
At2 = ktensor(A2);
% Gt = ktensor(G);

% err = test_compute_WB(A,G);
% % 
% err = test_compute_Wb(A,G);
% % 
% [err, err2] = test_compute_AFtAF(A,F);
% 
% [err, err2] = test_compute_B(A,F,G);

% err = test_SRMultV(A,ktensor(F));

err = test_SRMultM(At,At2);

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

function [err] = test_compute_Wb(A,G)

    d = ndims(A);
    rA = ncomponents(A);
    rG = ncomponents(G);
    N = sqrt(size(A,1));
    
    Wb = cell(rA,rG,d);
    lambda = zeros(rA,rG);
    
    tic
    for kk = 1:d
        A_temp = reshape(A.U{kk},N,N,rA);
        G_temp = G.U{kk};
        for ii = 1:rA
            for jj = 1:rG
                lambda(ii,jj) = A.lambda(ii)*G.lambda(jj);
                Wb{ii,jj,kk} = lambda(ii,jj)*A_temp(:,:,ii)'*G_temp(:,jj);    
            end
        end
    end
    toc
    
    W = als_sys_weights(A,G);
    Wb2 = W{2};
    
    err = zeros(rA,rA,d);
    for kk = 1:d
        for ii = 1:rA
            for jj = 1:rG        
                err(ii,jj,kk) = max(max(Wb{ii,jj,kk}-Wb2((ii-1)*N+(1:N),jj,kk)));
            end
        end
    end
    
end

function [err, err2] = test_compute_AFtAF(A,F)
    
    d = ndims(A); 
    R = size(F{1},2);
    rA = ncomponents(A);
    N = sqrt(size(A,1));
    
    
    tic
    AFtAF = cell(rA,rA,d);
    AUtAU = cell(rA,rA,d);
    for ii = 1:rA
        for jj = 1:rA
            for kk = 1:d
                A_temp = reshape(A.U{kk},N,N,rA); 
                AUtAU{ii,jj,kk} = (A_temp(:,:,ii)*F{kk})'*A_temp(:,:,jj)*F{kk};
            end
        end
    end
    
    for ii = 1:rA
        for jj = 1:rA
            for kk = 1:d
                AFtAF{ii,jj,kk} = 1;
                for nn = [1:kk-1 kk+1:d] 
                    AFtAF{ii,jj,kk} = AFtAF{ii,jj,kk}*AUtAU{ii,jj,nn};
                end
            end
        end
    end
    toc
    
    tic
    A_U = cell(d,1);
    AU = cell(d,1);
    AUtAU2 = zeros(rA*R,rA*R,d);
    AFtAF2 = cell(d,1);
    for nn = 1:d
        A_U{nn} = reshape(A.U{nn},N,N*rA);
        AU{nn} = blockTransposeV2H(blockTransposeH2V(A_U{nn},N)*F{nn},N);
        AUtAU2(:,:,nn) = AU{nn}'*AU{nn}; 
    end
    for kk = 1:d
        AFtAF2{kk} = prod(AUtAU2(:,:,[1:kk-1 kk+1:d]),3);
    end
    toc
    
    err = zeros(rA,rA,d);
    for ii = 1:rA
        for jj = 1:rA 
            for kk = 1:d   
                max_val = max(max([AFtAF{ii,jj,kk} AFtAF2{kk}((ii-1)*R+(1:R),(jj-1)*R+(1:R))]));
                err(ii,jj,kk) = max(max(AFtAF{ii,jj,kk} - AFtAF2{kk}((ii-1)*R+(1:R),(jj-1)*R+(1:R))))/max_val;
            end
        end
    end
    
    err2 = zeros(rA,rA,d);
    for ii = 1:rA
        for jj = 1:rA 
            for kk = 1:d   
                max_val = max(max([AUtAU{ii,jj,kk} AUtAU2((ii-1)*R+(1:R),(jj-1)*R+(1:R),kk)]));
                err2(ii,jj,kk) = max(max(AUtAU{ii,jj,kk} - AUtAU2((ii-1)*R+(1:R),(jj-1)*R+(1:R),kk)))/max_val;
            end
        end
    end
end

function [err, err2] = test_compute_B(A,F,G)

    d = ndims(A);   
    R = size(F{1},2);
    rA = ncomponents(A);
    rG = ncomponents(G);
    N = sqrt(size(A,1));
    
    W = als_sys_weights(A,G);
    WB2 = W{1};
    Wb = W{2};
    
    A_U = cell(d,1);
    AU = cell(d,1);
    AUtAU2 = zeros(rA*R,rA*R,d);
    AUtG = zeros(rA*R,rG,d);
    AFtAF2 = cell(d,1);
    for nn = 1:d
        A_U{nn} = reshape(A.U{nn},N,N*rA);
        AU{nn} = blockTransposeV2H(blockTransposeH2V(A_U{nn},N)*F{nn},N);
        AUtAU2(:,:,nn) = AU{nn}'*AU{nn}; 
        AUtG(:,:,nn) = AU{nn}'*G.U{nn};
    end
    tic
    B = cell(R,R,d);
    b = cell(R,R,d);
    for kk = 1:d
        AFtAF2{kk} = prod(AUtAU2(:,:,[1:kk-1 kk+1:d]),3);
        Z = prod(AUtG(:,:,[1:kk-1 kk+1:d]),3);
        for ii = 1:R
            for jj = 1:R
                B{ii,jj,kk} = 0;
                b{ii,jj,kk} = 0;
                for iii = 1:rA
                    for jjj = 1:rA
                        B{ii,jj,kk} = B{ii,jj,kk} + WB2((iii-1)*N+(1:N),(jjj-1)*N+(1:N),kk)*AFtAF2{kk}((iii-1)*R+ii,(jjj-1)*R+jj);                       
                    end
                    for jjj = 1:rG
                        b{ii,jj,kk} = b{ii,jj,kk} + Wb((iii-1)*N+(1:N),jjj,kk)*Z((iii-1)*R+ii,jjj);
                    end
                end
            end
        end
    end
    toc
    
    tic
    B2 = cell(d,1);
    b2 = cell(d,1);
    Y = cell(d,1);
    for n = 1:d
        Y{n} = prod(AUtAU2(:,:,[1:n-1 n+1:d]),3);
        Z = prod(AUtG(:,:,[1:n-1 n+1:d]),3);
        B2{n} = 0;
        b2{n} = 0;
        for ii = 1:rA
            for jj = 1:rA
                B2{n} = B2{n} + kron(Y{n}((ii-1)*R+(1:R),(jj-1)*R+(1:R)),W{1}((ii-1)*N+(1:N),(jj-1)*N+(1:N),n));
            end
            utemp = W{2}((ii-1)*N+(1:N),:,n)*Z((ii-1)*R+(1:R),:)';
            b2{n} =  b2{n} + utemp(:);
        end
    end
    toc
    
    err = zeros(R,R,d);
    for ii = 1:R
        for jj = 1:R
            for kk = 1:d   
                max_val = max(max([B{ii,jj,kk} B2{kk}((ii-1)*N+(1:N),(jj-1)*N+(1:N))]));
                err(ii,jj,kk) = max(max(B{ii,jj,kk} - B2{kk}((ii-1)*N+(1:N),(jj-1)*N+(1:N))))/max_val;
            end
        end
    end
    
    err2 = zeros(R,R,d);
    for ii = 1:R
        for jj = 1:1
            for kk = 1:d   
                max_val = max(max([b{ii,jj,kk} b2{kk}((ii-1)*N+(1:N),jj)]));
                err2(ii,jj,kk) = max(max(b{ii,jj,kk} - b2{kk}((ii-1)*N+(1:N),jj)))/max_val;
            end
        end
    end
end

function [err] = test_SRMultV(A,F)

    tic
    rF = ncomponents(F);
    rA = ncomponents(A);
    nd = ndims(F);

    % tensor access is slow, so do this at the expense of mem
    Alambda = A.lambda;
    Flambda = F.lambda;
    AU = cell(nd,1);
    FU = cell(nd,1);
    for k = 1:nd
        n = size(F.U{k},1);
        AU{k} = reshape(A.U{k},n,n,rA);
        FU{k} = F.U{k};
    end


    first = 1;
    for kA = 1:rA
        for kF = 1:rF
            clear T U
            U = cell(1,nd);
            for k = 1:nd
                U{k} = AU{k}(:,:,kA)*FU{k}(:,kF);
            end
            T = ktensor(Alambda(kA)*Flambda(kF),U);
            if first
                res = T;
                first = 0;
            else
                res = res + T;
            end
        end
    end
    toc
    
    tic
    res2 = SRMultV(A,F);
    toc
    
    err = norm(res-res2)/norm(res);
end

function [err] = test_SRMultM(A,B)

    tic
    rB = ncomponents(B);
    rA = ncomponents(A);
    nd = ndims(A);

    Alambda = A.lambda;
    Blambda = B.lambda;

    for k = 1:nd
        nA = round(sqrt(length(A.U{k}(:,1))));
        nB = round(sqrt(length(B.U{k}(:,1))));

        AU{k} = reshape(A.U{k},nA,nA,rA);
        BU{k} = reshape(B.U{k},nB,nB,rB);
    end

    first = 1;
    for kA = 1:rA
        for kB = 1:rB
            clear T U mat
            U = cell(1,nd);
            for k = 1:nd
                mat = AU{k}(:,:,kA)*BU{k}(:,:,kB);
                U{k} = mat(:);
            end
            T = ktensor(Alambda(kA)*Blambda(kB),U);
            if first
                res = T;
                first = 0;
            else
                res = res + T;
            end
        end
    end
    toc
    
    tic
    res2 = SRMultM(A,B);
    toc
    
    err = norm(res-res2)/norm(res);
end