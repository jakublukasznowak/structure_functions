
function [S,U] = edr_fit(S,var,dr,fitting_range,options)


arguments
    S (:,1) struct
    var (1,1) string
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty}
    fitting_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty}
    options.Scaling (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 2/3
    options.Factor (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 2.0
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'logmean'
    options.FittingPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 10
    options.Uncertainty (:,1) struct = []
end

if ~isempty(options.Uncertainty)
    U = options.Uncertainty;
else
    U = struct([]);
end

Nseg = size(S,1);

for i = 1:Nseg
    
    S(i).(strcat('fit_',var)) = fitting_range;
    S(i).(strcat('slp_',var)) = options.Scaling;
    S(i).(strcat('con_',var)) = options.Factor;

end






for i = 1:Nseg
    
    fit_range = S(i).(strcat('fit_',var));
    slp       = S(i).(strcat('slp_',var));
    C         = S(i).(strcat('con_',var));

    
    
    r = (1:length(S(i).(var)))*dr;

    ind_range = find( r>=fit_range(1) & r<=fit_range(2) );
    

    if strcmp(options.Method,'logmean')
        [rv_fit,sfc_fit] = logmean(r(ind_range),S(i).(var)(ind_range),options.FittingPoints);
    else
        rv_fit = r(ind_range);
        sfc_fit = S(i).(var)(ind_range);
    end

    
    % Fit (1): fixed slope
    
    offsetFixed = mean(log(sfc_fit)-slp*log(rv_fit));
    e_offsetFixed = std(log(sfc_fit)-slp*log(rv_fit))/sqrt(length(rv_fit));
    edrFixed = (exp(offsetFixed)/C)^(1/slp);

    S(i).(strcat('edr_',var)) = edrFixed;
    
    
    % Fit (2): free slope
    
    [p,s] = polyfit(log(rv_fit),log(sfc_fit),1);
    covM = (inv(s.R)*inv(s.R)')*s.normr^2/s.df;

    slpFree = p(1);
    e_slpFree = sqrt(covM(1,1));
    offsetFree = p(2);
    e_offsetFree = sqrt(covM(2,2));
    edrFree = (exp(offsetFree)/C)^(1/slpFree);

    S(i).(strcat('edr_',var,'_free')) = edrFree;
    S(i).(strcat('slp_',var,'_free')) = slpFree;
    
    
    % Linear correlation
    
    corM = corrcoef(log(rv_fit),log(sfc_fit));
    S(i).(strcat('r2_',var)) = corM(1,2);
    
        
    % Uncertainty
    
    U(i).(strcat('edr_',var)) = edrFixed/slp*e_offsetFixed;
    U(i).(strcat('edr_',var,'_free')) = edrFree/slpFree * sqrt( e_offsetFree^2 + (e_slpFree*log(edrFree))^2 );
    U(i).(strcat('slp_',var,'_free')) = e_slpFree;
    
end


% Combined 3D EDR
if all(ismember({'edr_ww','edr_uu','edr_vv'},fieldnames(S)))
    for i=1:Nseg
        S(i).edr_s2 = mean( [S(i).edr_ww, S(i).edr_uu, S(i).edr_vv] );
    end
end
if all(ismember({'edr_ww','edr_uu','edr_vv'},fieldnames(U)))
    for i=1:Nseg
        U(i).edr_s2 = mean( [U(i).edr_ww, U(i).edr_uu, U(i).edr_vv] );
    end
end


end