function [op,conv,diff] = makeop(fTens,BTens,noise_covTens,qTens,D,D2,lambda)
%MAKEOP creates the operator of the linear HJB equation.
% (ktensor cell := cell array with ktensors as elements.)
% Input:
%   fTens - f as a ktensor cell.
%   BTens - B as a ktensor cell.
%   noise_covTens - noise_cov as a matrix ktensor.
%   qTens - q as ktensor.
%   D - cell with first order differentiation operators. D{i} is
%   differentiation with respect to dimension i.
%   D2 - cell with second order differentiation operators. D2{i} is
%   differentiation with respect to dimension i.
%   lambda - the scalar lambda.
%Output:
%   op - the operator as a ktensor.
%   conv - vector corresponding to which dimensions that have a convection
%   term. conv(i) = 1 if dimension i has convection term, otherwise 0.
%   diff - vector corresponding to which dimensions that have a diffusion
%   term. diff(i) = 1 if dimension i has diffusion term, otherwise 0.
%
% See also MAIN_RUN.

% Elis Stefansson, Aug 5 2015

% the operator is divided into three terms:
% 1. left operator : -1/lambda*q
% 2. middle operator : f'*(gradOp)
% 3. right operator : 0.5*Tr(Hess*Sigma_t)

sizef = size(fTens);
d = sizef(1);

%% create left operator
opLeft = -1/lambda*DiagKTensor(qTens{1});

% debugging %
%save('opLeft','opLeft')
% debugging %

% debugging %
%opLeft_tol = opLeft;
%load('InvPenPerOrg_opLeft')
%opLeftdiff = norm(opLeft-opLeft_tol)
%norm(opLeft)
%error('done')
% debugging %

% debugging VTOL %
%fprintf('opleft')
%opLeft_tol = opLeft;
%load('opLeft_VTOL_org','opLeft')
%opLeft_org = opLeft;
%norm(opLeft_tol-opLeft_org)
%norm(opLeft_org-opLeft_org)
%norm(opLeft_org)
%error('done')
% debugging VTOL %

%% create middle operator
conv = ones(1,d); %start with ones.
opMid = [];
for i=1:d
    if norm(fTens{i}) == 0
        conv(i) = 0; %no convenction term in dimension i.
    else
        if isempty(opMid) == 1
            opMid = SRMultM(DiagKTensor(fTens{i}),D{i});
        else
            opMid = opMid+SRMultM(DiagKTensor(fTens{i}),D{i});
        end
    end
end

% debugging %
%save('opMid','opMid');
% debugging %

% debugging VTOL %
%fprintf('opMid')
%opMid_tol = opMid;
%load('opMid_VTOL_org','opMid')
%opMid_org = opMid;
%norm(opMid_tol-opMid_org)
%norm(opMid_org-opMid_org)
%norm(opMid_org)
%ncomponents(opMid_tol)
%ncomponents(opMid_org)
%error('done')
% debugging VTOL %

%% create right operator
diff = ones(1,d);

% create Hessian matrix
Hess = hessian(D,D2);

% debugging for VTOL %
%Hess_tol = Hess;
%load('hess_VTOL_org','hess')
%Hess_org = hess;
%sizeHess = size(Hess)
%for i = 1:sizeHess(1)
%    for j = 1:sizeHess(2)
%        [i,j]
%        norm(Hess_org{i,j}-Hess_tol{i,j})
%        norm(Hess_org{i,j}-Hess_org{i,j})
%        norm(Hess_org{i,j})
%    end
%end
%error('done')
% debugging for VTOL %

noise_covTensOp = DiagMatKTensor(noise_covTens);
BTensOp = DiagMatKTensor(BTens);
Sigma_t = TMatM(noise_covTensOp,BTensOp');
Sigma_t = TMatM(BTensOp,Sigma_t);

% debugging %
%Sigma_t_tol = Sigma_t;
%load('Sigma_t_VTOL_org','Sigma_t')
%Sigma_t_org = Sigma_t;
%sizeSigma = size(Sigma_t_org);
%for i = 1:sizeSigma(1)
%    for j = 1:sizeSigma(2)
%        [i,j]
%        norm(Sigma_t_org{i,j}-Sigma_t_tol{i,j})
%        norm(Sigma_t_org{i,j}-Sigma_t_org{i,j})
%        norm(Sigma_t_org{i,j})
%    end
%end
%error('done')
% debugging %

% create opRight = 0.5*Tr(Hess*Sigma_t)
opRight = [];
for i=1:d
    for j=1:d
        
        if i == j
            if norm(Sigma_t{j,i}) == 0
                diff(i) = 0; %no diffusion term in dimension i.
            end
        end
        
        % create term to right operator
        if norm(Sigma_t{j,i}) ~= 0
            if isempty(opRight) == 1
                opRight = SRMultM(Sigma_t{j,i},Hess{i,j});
            else
                opRight = opRight+SRMultM(Sigma_t{j,i},Hess{i,j});
            end
        end
        
    end
end
if isempty(opRight) == 0
    opRight = 0.5*opRight;
end

% debugging %
%save('opRight','opRight')
% debugging %

% debugging VTOL %
%fprintf('opRight')
%opRight_tol = opRight;
%load('opRight_VTOL_org','opRight')
%opRight_org = opRight;
%norm(opRight_tol-opRight_org)
%norm(opRight_org-opRight_org)
%norm(opRight_org)
%ncomponents(opRight_tol)
%ncomponents(opRight_org)
%error('done')
% debugging VTOL %

% debugging %
%h1 = 2*pi/200;
%h2 = 22/200;
%scale = h1^2*h2^2;
%opRight_tol = scale*opRight; %to have the same as them
%load('InvPenPerOrg_opRight')
%opRightdiff = norm(opRight-opRight_tol)
%norm(opRight_tol-opRight_tol)
%norm(opRight)
%error('done')
% debugging %

%% combine operator terms
if isempty(opMid) == 1
    if isempty(opRight) == 1
        op = opLeft;
    else
        op = opLeft+opRight;
    end
else
    if isempty(opRight) == 1
        op = opLeft+opMid;
    else
        op = opLeft+opMid+opRight;
    end
end


