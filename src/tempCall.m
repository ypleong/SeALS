function [ f2 ] = tempCall ( fD, t, x )
    % Call the symbolic function handle with arguments x1, ... xN using the
    % vector in x
    xarray = num2cell(x);
    f2 = feval(fD,xarray{:});
    
end

