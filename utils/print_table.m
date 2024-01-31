
function print_table(MOM,vars,ifcount,prec)

if nargin<4 || isempty(prec)
    prec = 2;
end
if nargin<3 || isempty(ifcount)
    ifcount = false;
end

% MOM.vstop = MOM.top - MOM.alt;
MOM.length_km = MOM.length/1000;

G = sortrows(groupsummary(MOM,"level",["mean","std"],union(vars,{'alt'})),"mean_alt",'descend');
T = groupsummary(MOM,[],["mean","std"],union(vars,{'alt'})); T.level = "all";
G = [G;T];

for i_l = 1:size(G,1)
    fprintf(' %20s',G.level(i_l))
    if ifcount
        fprintf(' & %5d',G.GroupCount(i_l))
    end
    for i_v = 1:numel(vars)
        fprintf(' & %6.*f (%.*f)',prec,G{i_l,['mean_',vars{i_v}]},prec,G{i_l,['std_',vars{i_v}]})
    end
    fprintf(' \\\\ \n')
end


end