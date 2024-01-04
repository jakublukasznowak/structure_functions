
function mask = ind2mask (indlist,L)

if nargin<2
    L = max(indlist(:));
end

mask = false(L,1);

for i = 1:size(indlist,1)
    mask(indlist(i,1):indlist(i,2)) = true;
end

end