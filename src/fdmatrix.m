function [DD] = fdmatrix( n, der, ord, per )
% FDMATRIX Finite difference matrix for der:th order derivative with ord:th 
% order accuracy in one dimension.
% Inputs:
%   n - grid size
%   ord - order of accuracy (interior, edge). If periodic: edge order of
%   accuracy will equal interior accuracy.
%   der - order of derivative.
%   per - periodic (1 - yes, [] - no).
% Outputs:
%   DD - the finite difference matrix.
%
%   See also MAKEDIFFOP.

% Elis Stefansson, Aug 9 2015

if nargin == 3
    per = 0;
end

if length(ord) < 2
    ord = [ord ord];
end

if mod(ord(1)+der,2) == 1 %ord(1)+der is odd
    
    c = fdcoeffF(der,0,-(ord(1)+der-1)/2:(ord(1)+der-1)/2);
        
    DD = make_fdmatrix(n,der,ord,per,c);
    
        
else %ord(1)+der is even
    
    c = fdcoeffF(der,0,-((ord(1)+der)/2-1):(ord(1)+der)/2); %gives the correct order of acc.
    c = [0 c]; %treat it as a central difference (compare to odd case)
    ord(1) = ord(1)+1; %this is a fake order, so we can use make_fdmatrix.
    
    DD = make_fdmatrix(n,der,ord,per,c);
        
end

end

function [DD] = make_fdmatrix(n, der, ord, per, c) % help function

dc1 = zeros(1,n);
dc1(1:(ord(1)+der-1)/2+1) = c((ord(1)+der-1)/2+1:end);
dc2 = zeros(1,n);
dc2(1:(ord(1)+der-1)/2+1) = c((ord(1)+der-1)/2+1:-1:1);
    
if per == 1
        
    add_to_dc1 = dc2(2:(length(c)+1)/2);
    add_to_dc1 = fliplr(add_to_dc1);
    dc1(end-length(add_to_dc1)+1:end) = add_to_dc1;
    
    add_to_dc2 = dc1(2:(length(c)+1)/2);
    add_to_dc2 = fliplr(add_to_dc2);
    dc2(end-length(add_to_dc1)+1:end) = add_to_dc2;
        
end
    
DD = toeplitz(dc2,dc1);

if per ~= 1
    for ii = 1:(ord(1)+der-1)/2
        cc = fdcoeffF(der,ii,1:ord(2)+der);
        cc = [cc zeros(1,n-length(cc))];
        DD(ii,:) = cc;
        DD(n-ii+1,:) = ((-1)^der)*cc(end:-1:1); %since "f(-x)".           
    end
end

end
