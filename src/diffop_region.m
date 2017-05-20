function [ Dpunch] = diffop_region( x, D, region, ord, Dedg, der)
%REGIONDIFFOP Adds boundary conditions to derivative for the region.
% Inputs:
%   x - the underlying grid. If x is a vector, then it is assumed the grid
%   is uniform along each dimension
%   D - the derivative operator
%   dim - number of nodes along each dimension
%   region - the interval where the derivative does not apply.
%   ord - order of accuracy in derivatives
%   Dedg - edge derivatives
% Outputs:
%   Dpunch - derivative operator incooporating region.
%
%   See also MAKEDIFFOP.

% Elis Stefansson, Aug 9 2015
% This is a new version incorporated for the Toolbox, explaining the 
% different notation.

d = length(x);
opDU = cell(d,d); %derivative operator in each direction
                  %column is derivative direction, row is U{i}
                  
if length(ord(:)) == 1 %same order for every dimension
    ord_old = ord;
    for i=1:d
        ord(i,1) = ord_old;
    end
else
    ord_old = ord(der,:);
    for i=1:d
        ord(i,:) = ord_old;
    end
end
                  
for ii = 1:d

    DCanv = addBCD_each( x{ii}, D{ii}, region(ii,:), ord(ii,:), Dedg{ii}, der);
    
    count = 0;
    if d == 1
         opDU{ii,ii} = DCanv(:);
    else
        for jj = 1:d           
            if ii == jj
                opDU{jj,ii} = [repmat(D{ii}(:),1,d-1) DCanv(:)];
            else
                id = eye(size(D{jj}));
                idRest = addBCD_each( x{jj}, id, region(jj,:), ord(jj,:), id, der);                
                opDU{jj,ii} = [repmat(id(:),1,count) idRest(:) repmat(id(:)-idRest(:),1,d-1-count)];
                count = count + 1;
            end         
        end
    end
end

Dpunch = cell(d,1);

for i=1:d
    Dpunch{i} = ktensor(opDU(:,i));
    Dpunch{i} = als2(Dpunch{i},10^-8);
end

end

function [ Dpunch] = addBCD_each( x, D, region, ord, Dedg, der)
    n = length(x);
    d_ind = find( (x >= region(1)) & (x <= region(2))); %nodes interior
    
    Dfill = zeros(n); %canvas for new D which we will update
    
    left = min(d_ind) - 1; %find the boundaries of the region
    right = max(d_ind) + 1;
    
    if mod(ord(1)+der,2) == 0
        ord(1) = ord(1)+1; %treat it as odd, like in fdmatrix
    end
    
    % special case when no edge derivates for region are needed.
    % (OBS: ord+der must be odd and can't be 1)
    if ord(1)+der == 3
        
        Drem = D;
        Drem(d_ind,:) = 0; %knock out the parts affecting the interior boundary
        
        Dpunch = Dfill+Drem;
        
    else
        
        if length(ord(:)) == 1
            ord = [ord ord]; %like in fdmatrix
        end
        
        d_egd = [(left+1)*ones((ord(1)+der-3)/2,1)-((ord(1)+der-3)/2:-1:1)'; (right-1)*ones((ord(1)+der-3)/2,1)+(1:(ord(1)+der-3)/2)'];
        
        %For the current direction, set up edge derivatives
        Dfill( (left+1-(ord(1)+der-3)/2):left, (left+1-(ord(2)+der-1)):(left+1) ) = Dedg(end-(ord(1)+der-3)/2:end-1,end-(ord(2)+der-1):end);
        Dfill(right:right+(ord(1)+der-5)/2,right-1:right-1+(ord(2)+der-1)) = Dedg(2:2+(ord(1)+der-5)/2,1:1+(ord(2)+der-1));

        Drem = D;
        Drem([d_egd; d_ind],:) = 0; %knock out the parts affecting the interior boundary
        
        Dpunch = Dfill + Drem;
        
    end
end
