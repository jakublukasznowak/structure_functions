
function [Tv,thv,th] = Tvthvth (thl,qv,ql,p)

% Constants
Rd = 287.04;  % [J/kg/K]
Rv = 461.5;
ep = Rd/Rv;
cpd= 1005;    % [J/kg/K]
Lv = 2.5e6;   % [J/kg]
p_ref = 1000e2; % [Pa]


th = thl + (p_ref./p).^(Rd/cpd)*Lv/cpd.*ql;

thv = th.*(1+(1/ep-1)*qv-ql);

Tv = (p/p_ref).^(Rd/cpd).*thv;

end