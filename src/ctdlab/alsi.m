function [Bk,err,Anorm,ext] = alsi(A,rnk,tol,stucktol,delta,ext,varargin)
% Main ALS iteration.  See als.m for arguments

    D=ndims(A);
    
    % parse optional inputs and/or set defaults
    params=inputParser;
    params.addParamValue('B',[],@(x) (isequal(class(x),'ktensor') || ...
        isempty(x)))
    params.addParamValue('density',knnz(A)/knumel(A),...
        @(x) isscalar(x) && x>0);
    params.addParamValue('maxit',500,@(x) isscalar(x) && x>0);
    params.addParamValue('verbose',0,@isscalar);
    params.addParamValue('dimorder',1:D, @(x) isequal(sort(x),1:D));
    params.addParamValue('Anorm',[],@(x) (isscalar(x) && x>=0) || isempty(x));
    params.addParamValue('errtype','rel',@(x) (ischar(x) && (isequal(x,'rel') ...
        || isequal(x,'abs'))));
    params.parse(varargin{:});

    % copy params
    B=params.Results.B;
    density=params.Results.density;
    maxit=params.Results.maxit;
    verbose=params.Results.verbose;
    dimorder=params.Results.dimorder;
    Anorm=params.Results.Anorm;
    errtype=params.Results.errtype;
    
    % print parameters
    if verbose
        fprintf('***alsr.m: required parameters\n')
        fprintf('      tol = %g, stucktol = %g, delta = %g, rnk = %d\n\n',...
            tol,stucktol,delta,rnk)
        fprintf('***                     optional parameters:\n')
        disp(params.Results)
        fprintf('\n')        
    end
    
    % if using relative error and norm of A not supplied, compute it
    if isempty(Anorm) && isequal(errtype,'rel')
        Anorm = fnorm(A);
    end

    % initial guess
    if isempty(B)
        
        % if no initial guess provided, use a random sparse one
        B=cell(1,D);
        for d=1:D
            B{d} = sprandn(size(A,d),rnk,density);
        end
        
        % normalize the initial guess
        B = ktensor(1,B);
        B = normalize(B);
        B = B.U;
    else
        Bk = B;
        B = B.U;
        rnk = length(Bk.lambda);
    end
        
    % residual error vector
    err = zeros(1,maxit);
    
    % Main ALS iteration.
    % This implementation uses matlab vectorized operations as much as
    % possible.  See Kolda & Bader (2009) for vectorized formulation.
    for iter = 1:maxit
        
        for k = dimorder
         
            % list of "active" dimensions
            ia = dimorder;
            ia(find(ia==k))=[];
            
            % coefficient matrix to be inverted
            Y = ones(rnk,rnk);
            for d = ia
                Y = Y .* (B{d}'*B{d});
            end
            
            % regularize, alpha just larger than eps...
            Y = full(Y) + 1e-10 * eye(rnk);
            
            % compute pseudo-inverse
            [UY,DY,VY] = svd(full(Y));
            irng = find((diag(DY))>(DY(1,1)*1e-10)); % I'm not sure this is needed
            VY = VY(:,irng);
            UY = UY(:,irng);
            DY = DY(irng,irng);
            Ypinv = VY * diag(1./diag(DY)) * UY';
            
            % khatri-rao product.  This function is in the sandia tensor toolbox.
            Z = mttkrp(A,B,k);
            
            % multiply by right pseudoinverse
            Z = Z * Ypinv;

            % normalize
            lambda = sqrt(abs(sum(Z.^2,1)))';
            Z = Z * spdiags(1./lambda,0,rnk,rnk);

            % update
            B{k} = Z;
            Bk = ktensor(lambda,B);
            
            % extra info for analysis
            if ~isempty(ext)
              ext.ls_cond = [ext.ls_cond, max(abs(diag(DY))) / min(abs(diag(DY)))];
              ext.nsig_svals = [ext.nsig_svals, length(irng)];
              ext.lambda_l1 = [ext.lambda_l1, norm(lambda,1)];
              ext.rank_iter = [ext.rank_iter; [rnk iter k]];
              ext.Bnorm = [ext.Bnorm, full(norm(Bk))];
            end
            
        end
        

        
        % truncate small elements by sparsity factor
        Bk = trncel(Bk,delta);
        
        % check for convergence, print some info   
        err(iter) = norm(A-Bk);
        if isequal(errtype,'rel')
            err(iter) = err(iter) / Anorm;
        end
        
        % error "velocity", used to check if stuck
        if iter>1
            derr = abs(err(iter)-err(iter-1));
        else
            derr = [];
        end
        
        if verbose
          fprintf('rnk = %d, iter = %d, err = %g, derr = %g\n', ...
            rnk,iter,full(err(iter)),full(derr))
          if ~isempty(ext)
          fprintf('\nCond ............... = %e\n',ext.ls_cond(end))
          fprintf('Significant svals .... = %d\n', ext.nsig_svals(end));
          fprintf('L1 norm of sval vector = %e\n', ext.lambda_l1(end));
          fprintf('||B||_F = ............ = %e\n\n', ext.Bnorm(end));
          end
        end
        
        % if converged, exit
        if err(iter)<tol
            if verbose
                fprintf('Converged to rank %d\n', rnk)
            end
            err = err(1:iter);
            return
        end
        
        % if stuck, exit
        if iter>1
            if derr < stucktol
                if verbose
                    fprintf('Stuck!\n')
                end
                err = err(1:iter);
                return
            end
        end
    end
    
    if verbose
        fprintf('Reached maximum iterations.\n\n')
    end
    
    return
end