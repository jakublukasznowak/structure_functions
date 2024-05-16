
function indlist = mask2ind (logvec,nan_replace)

if nargin<2 || isempty(nan_replace)
    nan_replace = false;
end

% MASK2IND converts a mask from logical vector into a list of start/end
% indices.
%
% INDLIST = mask2ind (LOGVEC) computes a 2-column matrix of indices
% denoting the start (1st column) and the end (2nd column) of the event
% denoted as TRUE values in the input logical vector LOGVEC.
%
% If the LOGVEC starts with TRUE, INDLIST(1,1) is 1.
% If the LOGVEC ends with TRUE, INDLIST(2,end) is length(LOGVEC).


logvec(isnan(logvec)) = nan_replace;

logvec2 = [false; logvec; false];


if max(logvec2)>0
    ip = find(diff(logvec2)== 1); 
    im = find(diff(logvec2)==-1); 
    indlist = [ip im-1];
else
    indlist = [];
end

    
end