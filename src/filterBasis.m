function [ Tout ] = filterBasis( Tin, filterIndex  )
% Retrieves the same tensor with the specified basis by filterIndex.
% Use it to visualize the contribution of each basis. e.g.
%
% nLambda = length(Tin.lambda);
% for i=1:nLambda
%     figure(i)
%     filterVector = zeros(nLambda,1);
%     filterVector(i) = 1;
%     plotkTensor(filterBasis(Tin,filterVector),gridT);
%     saveas(gcf,['basis_',num2str(i),'.png'])
% end

% Tin: input tensor
% filterIndex: vector containing 0 or 1 to indicate which basis to keep
    
    validateattributes(Tin,{'ktensor'},{'nonempty'})    
    dim = ndims(Tin);
    Uout = cell(dim,1);
    nLambda = length(Tin.lambda);
    validateattributes(filterIndex,{'double'},{'numel',nLambda})    
    index = 0;
    for i=1:nLambda
        if (filterIndex(i))
            index = index +1;
            lambdaOut(index) = Tin.lambda(i);
            for k=1:dim
                Uout{k}(:,index) = Tin.U{k}(:,i);            
            end
        end
    end
    Tout = ktensor(lambdaOut',Uout);


end

