function [ op ] = create_FP_op( fFPE, q, dt, D, D2,gridT, tol_err_op, explicit,x)
% 
% fTens, vector ktensor 
% qTens, matrix ktensor
% D, from makediffop
% D2, from makediffop
% dt, time difference (scalar) ktensor
%
% Write the operator H = Idendity + dt=(f*d/dx_i + q/2*d^2/dx_i/dx_j)

%% 
    
    fFPEdiff = jacobian(fFPE,x);
    

    fTens  = fcell2ftens( fsym2fcell(sym(fFPE) ,x), gridT);
    f_pTens = fcell2ftens( fsym2fcell(diag(fFPEdiff) ,x), gridT);
    if ( explicit == 1)
        dtTens = fcell2ftens( fsym2fcell(sym(dt)     ,x), gridT);
    else 
        dtTens = fcell2ftens( fsym2fcell(sym(-dt)     ,x), gridT);
    end
    qTens  = fcell2ftens( fsym2fcell(sym(q)      ,x), gridT);
    dtTen = DiagKTensor(dtTens{1}); 
    
    sizef = size(fTens);
    dim = sizef(1);
    
    %% create (fi*p),i ::it assumes f_i constant so it is equivalennt to f_i*p,i. If f_i is not constant add f_i,i*p at the end
    %% TODO: add f_i,i for the general case.
    
    conv = ones(1,dim); %start with ones.
    fipi = [];
    for i=1:dim
        if norm(fTens{i}) == 0
            conv(i) = 0; %no convenction term in dimension i.
        else
            if isempty(fipi) == 1
                fipi = SRMultM(DiagKTensor(fTens{i}),D{i});
            else
                fipi = fipi+SRMultM(DiagKTensor(fTens{i}),D{i});
            end
        end
    end
    
    %% Create f_i,i*p
    fiip = [];
    for i=1:dim
        if norm(f_pTens{i}) == 0
            conv(i) = 0; %no convenction term in dimension i.
        else
            if isempty(fiip) == 1
                fiip = DiagKTensor(f_pTens{i});
            else
                fiip = fiip + DiagKTensor(f_pTens{i});
            end
        end
    end
    %% create (q_ij*p),ij ::it assumes q_ij constant so it is equivalennt to q_ij*p,ij
    diff = ones(1,dim);

    % create Hessian matrix
    Hess = hessian(D,D2);

    noise_covTensOp = DiagMatKTensor(qTens);

    qijpij = [];
    for i=1:dim
        for j=1:dim


            % create term to right operator
            if norm(noise_covTensOp{j,i}) ~= 0
                if isempty(qijpij) == 1
                    qijpij = SRMultM(noise_covTensOp{j,i},Hess{i,j});
                else
                    qijpij = qijpij+SRMultM(noise_covTensOp{j,i},Hess{i,j});
                end
            end

        end
    end
    if isempty(qijpij) == 0
        qijpij = 0.5*qijpij;
    end

    %% combine operator terms 
    aop =  -fipi;
    if isempty(fiip) == 0
       aop = aop -fiip; 
    end
    if isempty(qijpij) == 0
       aop = aop + qijpij; 
    end
    icell = {{{@(x1)ones(size(x1)),1}}};
    iTens = fcell2ftens(icell,gridT);
    iTen = DiagKTensor(iTens{1}); 
    op_uncomp = iTen + SRMultM(dtTen, aop);
    
    %% Compress Operator
    fprintf('Attempt to compress operator, rank(op)=%d\n', ncomponents(op_uncomp));
    rank_op_uncomp = ncomponents(op_uncomp);
    tic;
    start_compress = tic;

    [op, err_op, iter_op, enrich_op, t_step_op, cond_op, noreduce] = als2(op_uncomp,tol_err_op);

    compress_time = toc(start_compress);
    toc;
    rank_op_comp = ncomponents(op);
    fprintf('Number of components after compression, %d\n', ncomponents(op));

end
