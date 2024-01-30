
function [S,U] = calc_edr(S,U,var,options)

arguments
    S (:,1) struct
    U (:,1) struct
    var (1,1) string
    options.Slope (1,1) {mustBeFinite, mustBeNonempty, mustBeReal} = 2/3;
    options.Constant (1,1) {mustBeFinite, mustBeNonempty, mustBeReal} = 2.0;
    options.FitRange (1,2) {mustBeFinite, mustBeNonempty, mustBePositive} = [10 60]; % m
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'direct'
    options.FitPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 10
end


Nseg = size(S,1);

for i_s = 1:Nseg
     
    fit_range  = options.FitRange;
    slopeFixed = options.Slope;
    C          = options.Constant;
    
    S(i_s).(strcat('fit_',var)) = options.FitRange;
    S(i_s).(strcat('slp_',var)) = options.Slope;
    S(i_s).(strcat('con_',var)) = options.Constant;

    
    r = (1:length(S(i_s).(var)))*S(i_s).dr;
    ind_range = find( r>=options.FitRange(1) & r<=fit_range(2) );

    rv = r(ind_range);
    sfc = S(i_s).(var)(ind_range);


    % Average SFC in log-equal bins (for LOGMEAN method)

    if strcmp(options.Method,'logmean')
        [rv_fit,sfc_fit] = logmean(rv,sfc,options.FitPoints);
    else
        rv_fit = rv;
        sfc_fit = sfc;
    end


    % Fit (1): fixed slope

    offsetFixed = mean(log(sfc_fit)-slopeFixed*log(rv_fit));
    edrFixed = (exp(offsetFixed)/C)^(1/slopeFixed);

    e_offsetFixed = std(log(sfc_fit)-slopeFixed*log(rv_fit))/sqrt(length(rv_fit));
    e_edrFixed = edrFixed/slopeFixed*e_offsetFixed;

    S(i_s).(strcat('edr_',var)) = edrFixed;
    U(i_s).(strcat('edr_',var)) = e_edrFixed;


    % Fit (2): free slope

    [p,s] = polyfit(log(rv_fit),log(sfc_fit),1);

    slopeFree = p(1);
    offsetFree = p(2);
    edrFree = (exp(offsetFree)/C)^(1/slopeFree);

    covM = (inv(s.R)*inv(s.R)')*s.normr^2/s.df;
    e_slopeFree = sqrt(covM(1,1));
    e_offsetFree = sqrt(covM(2,2));
    e_edrFree = edrFree/slopeFree * sqrt( e_offsetFree^2 + (e_slopeFree*log(edrFree))^2 );

    S(i_s).(strcat('edr_',var,'_free')) = edrFree;
    S(i_s).(strcat('slp_',var,'_free')) = slopeFree;
    U(i_s).(strcat('edr_',var,'_free')) = e_edrFree;
    U(i_s).(strcat('slp_',var,'_free')) = e_slopeFree; 


    % Linear correlation

    corM = corrcoef(log(rv_fit),log(sfc_fit));
    S(i_s).(strcat('r2_',var)) = corM(1,2);
        
    
end


% Combined EDR

if all(ismember({'edr_ww','edr_uu','edr_vv'},fieldnames(S)))
    for i_s = 1:Nseg
        S(i_s).edr_s2 = mean( [S(i_s).edr_ww, S(i_s).edr_uu, S(i_s).edr_vv] );
    end
end

if all(ismember({'edr_ww','edr_uu','edr_vv'},fieldnames(U)))
    for i_s = 1:Nseg
        U(i_s).edr_s2 = mean( [U(i_s).edr_ww, U(i_s).edr_uu, U(i_s).edr_vv] );
    end
end


end