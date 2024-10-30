datapath = '/mnt/mdisk/open/BOMEX/simulation_202407/';
addpath(genpath(myprojectpath))


% Variables
vars2load = {'thl','ql','qt'};
vars2full = {'thl','qt','ql','cf','z'};
varsthermo = {'thl','th','thv','qt','ql','qv','p','cf','T','Tv','ev','evs','rh','rho'};


% Constants
Rd = 287.04;  % [J/kg/K]
Rv = 461.5;
ep = Rd/Rv;
cpd= 1005;    % [J/kg/K]
Lv = 2.5e6;   % [J/kg]
T0 = 273.15;  % [K]
g  = 9.81;    % [m/s2]
p_ref = 1000e2; % [Pa]

ql_thresh = 0;

% Approximation of Clausius-Clapeyron
es = @(t) 6.112*exp(17.67*t./(t+243.5));  % hPa


% Read dimensions
info = ncinfo([datapath,'thl.nc'],'thl');
dims = num2cell(info.Size);
[Nx,Ny,Nz,Nt] = deal(dims{:});

% Create target struct
profhalf = cell2struct(repmat({nan(Nz,Nt)},1,numel(varsthermo)),varsthermo,2);
proffull = profhalf;


% Iterate over timesteps (not to overload memory)
for it = 1:Nt 
    
    % Load vars2load at half-levels and average in horizontal
    
    h.z = ncread([datapath,'thl.nc'],'z');
    h.time = ncread([datapath,'thl.nc'],'time');
    
    for iv = 1:numel(vars2load)
        var = vars2load{iv};
        x = ncread([datapath,var,'.nc'],var,[1 1 1 it],[Nx Ny Nz 1]);
        
        % Compute cloud fraction
        if strcmp(var,'ql')
            h.cf = x > ql_thresh;
        end
        
        h.(var) = squeeze(mean( x ,[1 2]));
        clear x
    end
    h.cf = squeeze(mean( h.cf ,[1 2]));
    

    % Calculate vars2full at full-levels
    
    f.time = h.time;
    for iv = 1:numel(vars2full)
        var = vars2full{iv};
        f.(var) = nan(size(h.(var)));
        f.(var)(2:end) = 0.5 * ( h.(var)(1:end-1) + h.(var)(2:end) );
    end
    
    % First full-level is z=0, i.e. sea surface
    f.thl(1) = 299.1;   
    f.qt(1)  = 22.45e-3; 
    f.ql(1)  = 0;
    f.z(1)   = 0;
    f.cf(1)  = 0;
    

    % Specific humidity
    
    h.qv = h.qt - h.ql;
    f.qv = f.qt - f.ql;
    
    
    % Compute p, Tv, thv, th by integrating hydrostatic balance from
    % the surface
    
    [f.p,f.th,f.thv,f.Tv,h.p,h.th,h.thv,h.Tv] = deal(nan(size(h.thl)));
    
    % 1st full-level
    f.p(1) = 1015e2;
    [f.Tv(1),f.thv(1),f.th(1)] = Tvthvth(f.thl(1),f.qv(1),f.ql(1),f.p(1));
    
    % 1st half-level
    p_old = f.p(1); p_new = 1000e2;
    while abs(p_new-p_old)>0.1
        Tv_old = Tvthvth(h.thl(1),h.qv(1),h.ql(1),p_old);
        Tv_mid = 0.5*(f.Tv(1)+Tv_old);
        p_old = p_new;
        p_new = f.p(1) * exp( -g*(h.z(1)-f.z(1)) / (Rd*Tv_mid) );
    end
    h.p(1) = p_new;
    [h.Tv(1),h.thv(1),h.th(1)] = Tvthvth(h.thl(1),h.qv(1),h.ql(1),h.p(1));
    
    % Integrate further in leapfrog spirit
    for iz = 2:Nz
        f.p(iz) = f.p(iz-1) * exp( -g*(f.z(iz)-f.z(iz-1)) / (Rd*h.Tv(iz-1)) );
        [f.Tv(iz),f.thv(iz),f.th(iz)] = Tvthvth(f.thl(iz),f.qv(iz),f.ql(iz),f.p(iz));
        
        h.p(iz) = h.p(iz-1) * exp( -g*(h.z(iz)-h.z(iz-1)) / (Rd*f.Tv(iz-1)) );
        [h.Tv(iz),h.thv(iz),h.th(iz)] = Tvthvth(h.thl(iz),h.qv(iz),h.ql(iz),h.p(iz));
    end
    
    
    % RH, Tv, density
    
    h.T = h.th.*(h.p/p_ref).^(Rd/cpd);  
    h.ev = h.p.*h.qv./(ep+h.qv*(1-ep));
    h.evs = es(h.T-T0)*100;
    h.rh = h.ev./h.evs;
    h.rho = h.p/Rd./h.Tv;
    
    f.T = f.th.*(f.p/p_ref).^(Rd/cpd);  
    f.ev = f.p.*f.qv./(ep+f.qv*(1-ep)); 
    f.evs = es(f.T-T0)*100;
    f.rh = f.ev./f.evs;
    f.rho = f.p/Rd./f.Tv;
       
    
    % Store
    
    fields = fieldnames(h);
    for iv = 1:numel(fields)
        var = fields{iv};
        if ismember(var,{'time','z'})
            profhalf.(var) = h.(var);
            proffull.(var) = f.(var);
        else
            profhalf.(var)(:,it) = h.(var);
            proffull.(var)(:,it) = f.(var);
        end
    end
    
end


save('bomex_profiles.mat','-struct','profhalf')