
function TURB = calc_thermo(TURB,MOM)


% Constants
Mv = 18.015;     % [g/mol]
Md = 28.9647;
R  = 8.31432;    % [J/mol/K]
Rd = R/Md*1000;  % [J/kg/K]
Rv = R/Mv*1000;
ep = Rd/Rv;
cpd= 1003.6;     % [J/kg/K]
% Lv = 2.5008e6;   % [J/kg]
T0 = 273.15;     % [K]
g  = 9.81;       % [m/s2]

% Approximation of Clausius-Clapeyron
es = @(t) 6.112*exp(17.67*t./(t+243.5));  % hPa


Nseg = size(TURB,1);

for i_s = 1:Nseg  
    p = MOM.MEAN_P(i_s)*100; % Pa
    
    th = TURB(i_s).T_DET + MOM.MEAN_THETA(i_s) + T0; % K potential temperature
    T = th * (p/1000/100)^(Rd/cpd);                  % K temperature
    
    
    % Water vapor
    
    rv = (TURB(i_s).MR_DET + MOM.MEAN_MR(i_s))/1000; % kg/kg mixing ratio
    qv = rv; %rv./(1+rv);                            % kg/kg mass fraction = specific humidity
    ev = p*rv./(rv+ep);                              % Pa water vapor partial pressure
    evs = es(T-T0)*100;                              % Pa saturation water vapor pressure
    rh = ev./evs;                                    % relative humidity
    
%     Tv = T.*(1+(1/ep-1)*qv);             % K virtual temperature
%     thv = Tv.*(1000*100./p).^(Rd/cpd);   % K virtual potential temperature
    thv = th.*(1+(1/ep-1)*qv);           % K virtual potential temperature
    
    
    % Save
           
    TURB(i_s).T    = T;         % K
    TURB(i_s).TH   = th;        % K
    TURB(i_s).THV  = thv;       % K
    
    TURB(i_s).RH   = rh*100;    % %
    
    TURB(i_s).B  = detrend( g/mean(thv) *(thv-mean(thv)) );                % buoyancy
    TURB(i_s).Bt = detrend( g/mean(thv) *(th -mean(th)) );                 % pseudobuoyancy theta-part
    TURB(i_s).Bq = detrend( g/mean(thv) *(1/ep-1)*(th.*qv-mean(th.*qv)) ); % pseudobuoyancy theta&q-part
    
end

end