function checklambda(G,B,noise_cov,R,lambda)
% CHECKLAMBDA checks that lambda is correctly specified for the special 
% case B=G and scalar noise_cov, R. Gives error if not.
% Inputs:
%   G
%   B
%   noise_cov
%   R
%   lambda
%   See MAIN_PROGRAM for specification.
%
%   See also MAIN_RUN.

% Elis Stefansson, Aug 6 2015

if isequaln(B,G) && length(noise_cov) == 1 && length(R) == 1
    
    % then lambda can be checked
    lambda_check = noise_cov*R;
    if lambda ~= lambda_check
        error('lambda is not correct')
    end
    
end
        
end
            