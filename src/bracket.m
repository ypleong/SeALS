function left_indices=bracket(xdata,xval)
% left_indices=bracket(xdata,xval) 
% Inputs:
% xdata - a vector of ascending ordered numbers
% xval - a point in the range of xdata
%
% Outputs:
% left_indices - index in xdata immediately to the left of xval

ndata=numel(xdata);
[xdata, ind] = sort(xdata);

left_indices=zeros(size(xval));
% Case 1
left_indices(  xval<xdata(2)  )=min(ind(1),ndata-1);
% Case 2
left_indices( xdata(ndata-1)<=xval  )=min(1, ind(ndata-1));
% Case 3
for k=2:ndata-2
    left_indices((xdata(k)<=xval) & (xval<xdata(k+1)))=ind(k);
end
if any(left_indices==0)
    error('bracket: not all indices set!')
end