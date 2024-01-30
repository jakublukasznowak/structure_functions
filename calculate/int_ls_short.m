
function LS = int_ls_short (x,y,options)

arguments
    x (:,1) {mustBeReal, mustBeNonempty}
    y (:,1) {mustBeReal, mustBeNonempty} = x
    options.Method (1,1) string {mustBeMember(options.Method,{'e-decay','integration'})} = 'e-decay'
    options.MaxLag (1,1) {mustBeInteger, mustBeFinite} = length(x)-1
    options.SubtractMean (1,1) logical = false
end


if options.SubtractMean
    xp = x - mean(x,'omitnan');
    yp = y - mean(y,'omitnan');
else
    xp = x;
    yp = y;
end

Fxy = mean(xp.*yp,'omitnan');


if strcmp(options.Method,'e-decay')
    lim = exp(-1);
elseif strcmp(options.Method,'integration')
    lim = 0;
end

lag = 0; lag1 = 0; lag2 = 0; 
rho1 = 1; rho2 = 1;

Lx = length(x);

while rho2>lim && lag<options.MaxLag
    
    lag = lag + 1;
    rho = mean( xp(lag+1:Lx) .* yp(1:Lx-lag), 'omitnan' )/Fxy;
    
    if ~isnan(rho)
        lag1 = lag2;
        rho1 = rho2;
        lag2 = lag;
        rho2 = rho;
    end
        
end

LC = interp1( [rho1,rho2], [lag1,lag2], lim );


if strcmp(options.Method,'integration')
    R = nan(lag2,1);
    for lag = 0:lag2-1
        R(lag+1) =  mean( xp(lag+1:Lx) .* yp(1:Lx-lag), 'omitnan' )/Fxy;
    end
    LS = trapz( [(0:lag2-1)';LC],[R;0] );
else
    LS = LC;
end


end
