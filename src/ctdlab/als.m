function [B,err,Anorm,ext] = als(A,tol,stucktol,delta,r0,rmax,varargin)
% Rank reduction by standard ALS.
%
% Required arguments:
%   A = sparse ktensor to reduce
%   tol = desired relative accuracy
%   stucktol = tolerance for being "stuck"
%   delta = all elements below this value are set to zero.  If [], no
%     truncation is performed.
%   r0 = initial separation rank.
%   rmax = maximum separation rank of approximation.  If [], defaults to
%     rank of A.
%
% Optional arguments:  See inputParser for defaults.
%   maxit = maximum number of iterations.
%   B = initial guess.  If empty, generate random sparse ktensor of rank r0
%   density = density of random factors.  If empty, use same density as A.
%   dimorder = D-vector, the order in which dimensions are fit in als.
%     Defaults to 1:D.
%   Anorm = norm of A.  Can re-use if previously computed, otherwise it is
%     computed within the main ALS iteration.
%   errtype = type of error to use.  'rel' for relative or 'abs' for
%     absolute.  Default is 'rel'.
%   verbose = 1, print iteration info.  0, none.
%
% Output:
%   B = low-rank approximation of A
%   err = error per ALS iteration
%   Anorm = norm of A.  You might want this if it was expensive to compute.
%   ext = strcuture containing info about the calculation
    
    D = ndims(A);
    
    % parse optional inputs and/or set defaults
    params=inputParser;
    params.addParamValue('B',[],@(x) (isequal(class(x),'ktensor') || ...
        isempty(x)))
    params.addParamValue('density',knnz(A)/knumel(A),...
        @(x) (isscalar(x) && x>0));
    params.addParamValue('maxit',100,@(x) isscalar(x) && x>0);
    params.addParamValue('verbose',0,@isscalar);
    params.addParamValue('dimorder',1:D, @(x) isequal(sort(x),1:D));
    params.addParamValue('Anorm',[],@(x) (isscalar(x) && x>=0)||isempty(x));
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
    
    % if the input is rank 1, exit without doing anything.
    if length(A.lambda)==1
        if verbose
            disp('Tensor already has rank 1, no reduction possible.')
        end
        B = A;
        err = 0;
        return
    end
    
    % start with rank 1 by default
    if isempty(r0)
      r0 = 1;
    end
    
    % if rmax is empty or larger than rank(A), set to rank(A)
    if isempty(rmax) && isempty(B)
        rmax = length(A.lambda);
    end
    if isempty(rmax) && ~isempty(B)
        rmax = length(B.lambda);
    end
    
    % if no initial guess specified, generate a random sparse one
    if ~isempty(B)
        r0 = length(B.lambda);
    end
    
    % print parameters
    if verbose
        fprintf('***als.m: Required parameters\n')
        fprintf('      tol = %g, stucktol = %g, delta = %g, rmax = %d\n\n',...
            tol,stucktol,delta,rmax)
        fprintf('***                   Optional parameters\n')
        disp(params.Results)
        fprintf('\n')
    end  
    
    % extra info for analysis
    ext.Bnorm = [];
    ext.ls_cond = [];    
    ext.lambda_l1 = [];
    ext.nsig_svals = [];
    ext.rank_iter = [];
    
    % Main iteration
    err = zeros(1,rmax-r0+1);
    for rnk = r0:rmax
        
        % if you haven't converged and rnk equals rank of A, exit.
        if (r0<rmax) && (rnk==length(A.lambda))
            % Message added by MR
            fprintf('Rank of reduction = rank of original matrix.\n')
            fprintf('Returning original matrix!\n')
            B = A;
            return
        end
        
        % track total iteration count (includes subcalls)
        iter = 1+rnk-r0;
        
        % if not the first iteration, increase the rank by 1 using a random
        % term with specified density
        if rnk > r0
            B = rankup(B);
        end

        % fit using main ALS iteration
        [B,err_r,Anorm,ext] = alsi(A,rnk,tol,stucktol,delta,ext,...
            'density',density,'B',B,'maxit',maxit,'verbose',verbose,...
            'dimorder',dimorder,'Anorm',Anorm,'errtype',errtype);
        err(iter) = err_r(end);
        
        
        % check for convergence
        if err(iter)<tol
            err = err(1:iter);
            return
        end
        
        % if at max rank, exit
        if rnk == rmax
            if verbose
                fprintf('***als.m: No tensor of rank rmax')
                fprintf(' found for desired accuracy.\n')
            end
            
            return
        end
        
    end

end

function A=rankup(A)
% Increase rank of spktensor by 1, using random component

    D=ndims(A);
    
    B=cell(1,D);
    for d=1:D
        B{d}=randn(size(A.U{d},1),1);
    end
    
    % Normalize the contributing term
    B = ktensor(1,B);
    B = normalize(B);
    
    % Set the norm to be similar to the tensor A
    % Note, I made the guess an order of magnitude smaller than the
    % smallest s-values because we hope for a fast decay in s-values -MR 
    B.lambda = min(abs(A.lambda))/10;
    
    A = A+B;
    A = arrange(A);
end