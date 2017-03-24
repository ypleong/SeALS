function [Flambda, FU, status] = als_onestep(Glambda,GU,Flambda,FU,alpha)
% ALS_ONESTEP performs one step of the ALS algorithm. Based on the version
% in Baylkin and Mohlenkamp 2005.
% Inputs:
%   F is the current iterate where:
%      Flambda - normailization constants for F.
%      FU - vectors (basis functions) for F.
%   G is the tensor begin approximated where:
%      Glambda - normailization constants for G.
%      GU - vectors (basis functions) for G.
%   alpha - the regularization constant.
% Outputs:
%   Updated F in form of:
%      Flambda (as before).
%      FU (as before).
%   status - 0 if ill-conditioned matrices. Otherwise 1.
%
% See also ALS.

%%% Old comments %%%
% merged to match als_onestep_a.m to avoid any confusion.
% this is now the primary branch.

% based on the version in Beylkin and Mohlenkamp 2005
% performs one step of the ALS algorithm


% F is the current iterate
% Glambda and GU represent the tensor being approximated
% alpha is the regularization
%%% Old comments %%%

status = 1;

nd = length(FU);
rG = length(Glambda);
rF = length(Flambda);


% tensor access is slow so this function now takes the inputs and outputs
% as lambda and U.

 % compute inner products needed for 3.3, really only need half of these
 % because the matrix is symmetric (in the real case)
 Fdot = zeros(rF,rF,nd);
 for n = 1:nd
	 Fdot(:,:,n) = FU{n}' * FU{n};
 end
 
 % compute inner products needed for 3.4
 GFdot = zeros(rG,rF,nd);
 for n = 1:nd
	 GFdot(:,:,n) = GU{n}' * FU{n};
 end

for k = 1:nd
	
    % update inner products needed for 3.3 and 3.4
	if k > 1
		Fdot(:,:,k-1) = FU{k-1}' * FU{k-1};
		GFdot(:,:,k-1) = GU{k-1}' * FU{k-1};
	end

	% construct the matrix B in (3.3)
	% the .' is to deal with the complex case, this has not been tested and
	% could be wrong. In the real case the matrix is symmetric and this will not
	% be an issue.
    idx = [1:(k-1) (k+1):nd];
	B = prod(Fdot(:,:,idx),3).';
	B = B + alpha*eye(size(B));
	
	% form product of G and F inner products
	GF = prod(GFdot(:,:,idx),3);
	
	b = (GU{k}*diag(Glambda)*GF).';
    
    if cond(B) > 1e13
        status = 0;
    end
    
    c = B\b;

    sF = sqrt(sum(c.^2,2));
    Flambda = sF;
    FU{k} = c.'*(diag(1./sF));
	
end

%%%%%%%%%%%%
% OLD CODE %
%%%%%%%%%%%%


% % based on the version in Beylkin and Mohlenkamp 2005
% % performs one step of the ALS algorithm
% % something may be wrong when the tensor F is of rank > 1 because in this case
% % the error can grow, which should not be possible. (This definitely seems to be
% % the problem...
% 
% 
% % F is the current iterate
% % G is the tensor being approximated 
% % alpha is the regularization
% G = arrange(G);
% nd = ndims(G);
% sizeG = size(G);
% rG = ncomponents(G);
% rF = ncomponents(F);
% 
% for k = 1:nd
% 	
%     % compute inner products needed for 3.3, really only need half of these
%     % because the matrix is symmetric (in the real case)
%     for n = 1:nd
%         Fdot(:,:,n) = F.U{n}' * F.U{n};
%     end
% 
%     % compute inner products needed for 3.4
%     for n = 1:nd
%         GFdot(:,:,n) = G.U{n}' * F.U{n};
%     end
% 
% 	% construct the matrix B in (3.3)
% 	% the .' is to deal with the complex case, this has not been tested and
% 	% could be wrong. In the real case the matrix is symmetric and this will not
% 	% be an issue.
% 	B = prod(Fdot(:,:,setdiff(1:nd,k)),3).';
% 	B = B + alpha*eye(size(B));
% 	
% 	% form product of G and F inner products
% 	GF = prod(GFdot(:,:,setdiff(1:nd,k)),3);
% 	
% 	% form vectors for 3.4 and solve the linear system in 3.5
% 	% this, or the update to F is probably where there is an issue
% 	
% 	for jk = 1:sizeG(k)
% 		b(:,jk) = zeros(rF,1);
% 		for lp = 1:rF
% 			for l = 1:rG
% 				b(lp,jk) = b(lp,jk) + G.lambda(l)*G.U{k}(jk,l)*GF(l,lp);
% 			end
% 		end
% 		c(:,jk) = B\b(:,jk);
% 	end
% 	
% 	for l = 1:rF
% 		sF = norm(c(l,:));
% 		F.lambda(l) = sF;
% 		for jk = 1:sizeG(k)
% 			F.U{k}(jk,l) = c(l,jk)/sF;
% 		end
%     end
% 	
%     F = arrange(F);
% end
