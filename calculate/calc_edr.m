
function [S,U] = calc_edr(S,U,vars,options)

arguments
    S (:,1) struct
    U (:,1) struct
    vars (1,:) cell
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'direct'
    options.FittingPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 10
end


Nseg = size(S,1);
Nvar = numel(vars);

for i_s = 1:Nseg
    
    for i_v = 1:Nvar
        
        var = vars{i_v};
        
        fit_range = S(i_s).(strcat('fit_',var));
        slp       = S(i_s).(strcat('slp_',var));
        C         = S(i_s).(strcat('con_',var));

        
        r = (1:length(S(i_s).(var)))*S(i_s).dr;

        ind_range = find( r>=fit_range(1) & r<=fit_range(2) );


        if strcmp(options.Method,'logmean')
            [rv_fit,sfc_fit] = logmean(r(ind_range),S(i_s).(var)(ind_range),options.FittingPoints);
        else
            rv_fit = r(ind_range);
            sfc_fit = S(i_s).(var)(ind_range);
        end


        % Fit (1): fixed slope

        offsetFixed = mean(log(sfc_fit)-slp*log(rv_fit));
        e_offsetFixed = std(log(sfc_fit)-slp*log(rv_fit))/sqrt(length(rv_fit));
        edrFixed = (exp(offsetFixed)/C)^(1/slp);

        S(i_s).(strcat('edr_',var)) = edrFixed;


        % Fit (2): free slope

        [p,s] = polyfit(log(rv_fit),log(sfc_fit),1);
        covM = (inv(s.R)*inv(s.R)')*s.normr^2/s.df;

        slpFree = p(1);
        e_slpFree = sqrt(covM(1,1));
        offsetFree = p(2);
        e_offsetFree = sqrt(covM(2,2));
        edrFree = (exp(offsetFree)/C)^(1/slpFree);

        S(i_s).(strcat('edr_',var,'_free')) = edrFree;
        S(i_s).(strcat('slp_',var,'_free')) = slpFree;


        % Linear correlation

        corM = corrcoef(log(rv_fit),log(sfc_fit));
        S(i_s).(strcat('r2_',var)) = corM(1,2);


        % Uncertainty

        U(i_s).(strcat('edr_',var)) = edrFixed/slp*e_offsetFixed;
        U(i_s).(strcat('edr_',var,'_free')) = edrFree/slpFree * sqrt( e_offsetFree^2 + (e_slpFree*log(edrFree))^2 );
        U(i_s).(strcat('slp_',var,'_free')) = e_slpFree;  
        
    end
    
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