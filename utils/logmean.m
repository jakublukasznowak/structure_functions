
function [xnew,ynew] = logmean (x,y,Nbin)

% LOGMEAN averages the values from two input vectors X and Y in NBIN bins
% which are log-equally distributed across the range of X values.

bin_edge = exp( linspace(log(min(x)),log(max(x)),Nbin+1)' );
bin_edge(1) = min(x);
bin_edge(end) = max(x);
bin_ind = discretize(x,bin_edge);

xnew = accumarray(bin_ind(:),x(:),[],@mean);
ynew = accumarray(bin_ind(:),y(:),[],@mean);

emptyBins = ~ismember((1:Nbin)',bin_ind);
Neb = sum(emptyBins);
if Neb > 0
    warning('LOGMEAN:EmptyBins','%d bins are empty.',Neb)
    xnew = xnew(~emptyBins);
    ynew = ynew(~emptyBins);
end

end