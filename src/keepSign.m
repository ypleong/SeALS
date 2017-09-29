function [ Tcell, flipSign] = keepSign( Tcell )
% keep the signs of the basis 

    nCell = length(Tcell);
    dim = ndims(Tcell{1});
    flipSign = zeros(nCell,1);
    for j=2:nCell
        for k=1:min(length(Tcell{j-1}.lambda),length(Tcell{j}.lambda))
            for i=1:dim
                uproduct = Tcell{j}.U{i}(:,k)'*Tcell{j-1}.U{i}(:,k);
                if ( uproduct < 0   )
                    flipSign(j) = flipSign(j) + 1;
                    Tcell{j}.U{i}(:,k)   = -Tcell{j}.U{i}(:,k);
                    %outTens{j}.U{i+1}(:,k) = -outTens{j}.U{i+1}(:,k);
                end                    
            end
            if mod(flipSign,2)~=0
                fprintf('Trouble fixing signs for mode %d\n', k);
            end
        end
    end
end
