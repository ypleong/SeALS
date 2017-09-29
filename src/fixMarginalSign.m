function [ outTens, sumK, sumIndexK ] = fixMarginalSign( iTens, gridT )
% fix the sign of each tensor to solve the pair of sign indeterminacy

    eps = 1e-7;

    dim = ndims(iTens);
    nLambda = length(iTens.lambda);
    outTens = iTens;
    for k=1:nLambda
        for i=1:2:dim
            sumIndex=1;
            sumi = mean(iTens.U{i}(1:round(end/sumIndex),k))/max(abs(iTens.U{i}(1:round(end/sumIndex),k)));
            while ( abs(sumi) < eps && sumIndex < 20 )
               sumIndex = sumIndex+1;
               sumi = mean(iTens.U{i}(1:round(end/sumIndex),k))/max(abs(iTens.U{i}(1:round(end/sumIndex),k)));
            end
            sumK(k) = sumi;
            sumIndexK(k) = sumIndex;
            if ( sumi <0  )
                outTens.U{i}(:,k)   = -iTens.U{i}(:,k);
                outTens.U{i+1}(:,k) = -iTens.U{i+1}(:,k);
            end                    

        end
    end
end

