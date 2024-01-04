
function indlist = mask2ind (logvec)

% MASK2IND converts a mask from logical vector into a list of start/end
% indices.
%
% INDLIST = mask2ind (LOGVEC) computes a 2-column matrix of indices
% denoting the start (1st column) and the end (2nd column) of the event
% denoted as TRUE values in the input logical vector LOGVEC.
%
% If the LOGVEC starts with TRUE, INDLIST(1,1) is 1.
% If the LOGVEC ends with TRUE, INDLIST(2,end) is length(LOGVEC).


logvec(isnan(logvec)) = 0;
Lc = length(logvec);

if max(logvec)>0
    ip = find(diff(logvec)== 1); Lp=length(ip);
    im = find(diff(logvec)==-1); Lm=length(im);
    if ip(1)<im(1)
        if Lp==Lm
            indlist = [ip+1 im];
        elseif Lp>Lm
            indlist = [ip+1 [im;Lc]];
        end
    else
        if Lp==Lm
            indlist = [[1;ip+1] [im;Lc]];
        elseif Lp<Lm
            indlist = [[1;ip+1] im];
        end
    end
else
    indlist=[];
end

    
end