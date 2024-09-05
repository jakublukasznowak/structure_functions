
function [S,U] = calc_edr(S,U,var,fit_range,C,options)

arguments
    S (:,1) struct
    U (:,1) struct
    var (1,1) string
    fit_range (1,2) {mustBeFinite, mustBeNonempty, mustBePositive}
    C (1,1) {mustBeFinite, mustBeNonempty, mustBeReal}
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'direct'
    options.FitPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 10
    options.Slope (1,1) {mustBeFinite, mustBeNonempty, mustBeReal} = 2/3;
    options.dC (1,1) {mustBeFinite, mustBeNonempty, mustBeReal} = 0;
end


slpFixed = options.Slope;


Nseg = size(S,1);

for i_s = 1:Nseg
    
    S(i_s).(strcat('fit_',var)) = fit_range;
    S(i_s).(strcat('con_',var)) = C;
    S(i_s).(strcat('slp_',var)) = options.Slope;
    
    
    % Select fit range
    
    r = (1:length(S(i_s).(var)))*S(i_s).dr;
    ind_r = find( r>=fit_range(1) & r<=fit_range(2) );
    rv = r(ind_r);
    sfc = S(i_s).(var)(ind_r);


    % Average SFC in log-equal bins (for LOGMEAN method)

    if strcmp(options.Method,'logmean')
        [rv_fit,sfc_fit] = logmean(rv,sfc,options.FitPoints);
    else
        rv_fit = rv;
        sfc_fit = sfc;
    end


    % Fit (1): fixed slope

    logO = mean(log(sfc_fit)-slpFixed*log(rv_fit));
    e_logO = std(log(sfc_fit)-slpFixed*log(rv_fit))/sqrt(length(rv_fit));
    
    edr = (exp(logO)/C)^(1/slpFixed);
    e_edr = sqrt( (edr/slpFixed*e_logO).^2 + (edr/slpFixed/C*options.dC).^2 );

    S(i_s).(strcat('edr_',var)) = edr;
    U(i_s).(strcat('edr_',var)) = e_edr;


    % Fit (2): free slope

    [p,s] = polyfit(log(rv_fit),log(sfc_fit),1);

    slp = p(1);
    logOstar = p(2);
    
    covarM = (inv(s.R)*inv(s.R)')*s.normr^2/s.df;
    e_slp = sqrt(covarM(1,1));
    e_logOstar = sqrt(covarM(2,2));
    
    S(i_s).(strcat('slp_',var,'_free')) = slp;
    U(i_s).(strcat('slp_',var,'_free')) = e_slp; 
    
    edrFree = (exp(logOstar)/C)^(1/slp);
    e_edrFree = edrFree/slp * sqrt( e_logOstar^2 + (e_slp*log(edrFree))^2 );
    
    S(i_s).(strcat('edr_',var,'_free')) = edrFree;
    U(i_s).(strcat('edr_',var,'_free')) = e_edrFree;
    

    % Linear correlation

    corM = corrcoef(log(rv_fit),log(sfc_fit));
    S(i_s).(strcat('R2_',var)) = corM(1,2);
        
    
end


% Combined EDR

if all(ismember({'edr_ww','edr_uu','edr_vv'},fieldnames(S)))
    for i_s = 1:Nseg
        S(i_s).edr_S2 = mean( [S(i_s).edr_ww, S(i_s).edr_uu, S(i_s).edr_vv] );
    end

    if all(ismember({'edr_ww','edr_uu','edr_vv'},fieldnames(U)))
        for i_s = 1:Nseg
            x = [S(i_s).edr_ww, S(i_s).edr_uu, S(i_s).edr_vv];
            s = [U(i_s).edr_ww, U(i_s).edr_uu, U(i_s).edr_vv];
%             xw = sum(x./s.^2)/sum(1./s.^2);
            xw = mean(x);
            S(i_s).edr_S2 = xw;
            U(i_s).edr_S2 = sqrt( max([1/sum(1./s.^2) sum( (x-xw).^2./s.^2 ) /2 / sum(1./s.^2)]) );
%             U(i_s).edr_S2 = mean( [U(i_s).edr_ww, U(i_s).edr_uu, U(i_s).edr_vv] );
        end
    end


end