
function TURB = calculate_thermo(TURB,MOM)

Mv = 18.015;     % [g/mol]
Md = 28.9647;
R  = 8.31432;    % [J/mol/K]
Rd = R/Md*1000;  % [J/kg/K]
Rv = R/Mv*1000;
ep = Rd/Rv;
cpd= 1003.6;     % [J/kg/K]
T0 = 273.15;     % [K]
g  = 9.81;       % [m/s2]

% Approximation of Clausius-Clapeyron
es = @(t) 6.112*exp(17.67*t./(t+243.5));  % hPa


Nseg = size(TURB,1);

for i_s = 1:Nseg  
    p = MOM.MEAN_P(i_s)*100; % Pa
    th = TURB(i_s).T_DET + MOM.MEAN_THETA(i_s) + T0; % K potential temperature
    q = (TURB(i_s).MR_DET + MOM.MEAN_MR(i_s))/1000; % kg/kg mixing ratio
        
    T = th * (p/1000/100)^(Rd/cpd); % K temperature
    thv = th.*(1+(1/ep-1)*q); % K virtual potential temperature
    
    ev = p*q./(q+ep); % Pa water vapor partial pressure
    evs = es(T-T0)*100; % Pa saturation water vapor pressure
    rh = ev./evs; % relative humidity
           
    
    TURB(i_s).T    = T;
    TURB(i_s).TH   = th;
    TURB(i_s).TH_V = thv; 
    TURB(i_s).RH   = rh*100;
    
    TURB(i_s).B  = detrend( g/mean(thv)*(thv-mean(thv)) ); % buoyancy
    TURB(i_s).Bt = detrend( g/mean(th) *(th -mean(th)) );  % pseudobuoyancy t-part
    TURB(i_s).Bq = detrend( g/mean(th) *(1/ep-1)*(th.*q-mean(th.*q)) ); % pseudobuoyancy tq-part
end

end