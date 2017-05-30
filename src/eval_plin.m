function yval = eval_plin ( xdata, ydata, xval )
% yval = eval_plin ( xdata, ydata, xval )
%
% Inputs:
% xdata - a n-length vector of ascending ordered numbers 
% ydata - values at xdata (n x r) r is the number of columns
% xval - a point in the range of xdata
% 
% Outputs:
% yval - interpolated value at xval for each column

left_index = bracket(xdata,xval);

if left_index ~= 1 && left_index ~= length(xdata)-1
    if (xdata(1) > xdata(end) && xdata(left_index) < xval) ...
            || (xdata(1) < xdata(end) && xdata(left_index) > xval)
        left_index = left_index-1;
    end   
end    
yval = ydata(left_index,:) + (xval - xdata(left_index,:))...
    *(ydata(left_index+1,:) - ydata(left_index,:))/(xdata(left_index+1) - xdata(left_index));

return