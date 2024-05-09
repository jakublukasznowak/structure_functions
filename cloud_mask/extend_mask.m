
function mask_new = extend_mask (mask,n,method)

% Extend each event marked by TRUE values in a logical MASK by its N width
% in front of and behind this event.


L = length(mask);

indlist = mask2ind(mask);

if strcmp(method,'constant') && size(indlist,1)>0
    
    indlist(:,1) = indlist(:,1) - n;
    indlist(:,2) = indlist(:,2) + n;
    
elseif strcmp(method,'proportional') && size(indlist,1)>0
    
    d = indlist(:,2) - indlist(:,1);
    indlist(:,1) = indlist(:,1) - n*d;
    indlist(:,2) = indlist(:,2) + n*d;
    
end

indlist(indlist<1) = 1;
indlist(indlist>L) = L;

mask_new = double(ind2mask(indlist,L));

mask_new(isnan(mask)) = nan;

end