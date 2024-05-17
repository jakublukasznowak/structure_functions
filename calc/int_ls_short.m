
% Calculate integral length scale LS [in points] based on the correlation
% function of data vectors X and Y (i.e. autocorrelation if Y=X,
% crosscorrelation otherwise).
%
% Optional Name-Value arguments:
% 'Method':
%       'e-decay'      -  (default) LS is the crossing of exp(-1) level
%       'zero'         -  LS is the crossing of the first zero
%       'integrate'    -  LS is the integral of correlation function up to
%                         the first zero crossing
%       'cum-integrate' - LS is the maximum of cumulative integral of the
%                         correlation function
% 'MaxLag': maximum lag to consider [in points]
% 'MaxLagFactor': max lag definition as factor x exp(-1) crossing for
%   integration methods
% 'dr': displacement between sample points [in meters]
% 'SubstractMean': whether to substract mean from X and Y before
%   calculating correlation function (true/false)
% 'Plot': selects whether to show a diagnostic plot (true/false)


function [LS,fig] = int_ls_short (x,y,options)

arguments
    x (:,1) {mustBeReal, mustBeNonempty}
    y (:,1) {mustBeReal, mustBeNonempty} = x
    options.Method (1,1) string {mustBeMember(options.Method,...
        {'e-decay','zero','integrate','cum-integrate'})} = 'e-decay'
    options.MaxLag (1,1) {mustBeInteger, mustBeFinite} = length(x)-1
    options.MaxLagFactor (1,1) {mustBePositive} = Inf
    options.dr (1,1) {mustBePositive} = 1
    options.SubtractMean (1,1) logical = false
    options.Plot (1,1) logical = false
end


% Substract mean
if options.SubtractMean
    xp = x - mean(x,'omitnan');
    yp = y - mean(y,'omitnan');
else
    xp = x;
    yp = y;
end

% (Co)variance
Fxy = mean(xp.*yp,'omitnan');

% Select critical limit value
if strcmp(options.Method,'e-decay')
    lim = exp(-1);
else
    lim = 0;
end


% Iterate calculating correlation function until the limit is reached

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


% Interpolate to get precise limit crossing position

if rho2>lim
    LC = options.MaxLag;
    
    if ismember(options.Method,{'e-decay','zero'})
        warning('INT_LS_SHORT:NoCrossing','\nNo %.1f crossing. Take max lag %.1f.',lim,LC)
    elseif strcmp(options.Method,'integrate')
        warning('INT_LS_SHORT:NoCrossing','\nNo %.1f crossing. Integrating to max lag %.1f.',lim,LC)
    end
    
else
    LC = interp1( [rho1,rho2], [lag1,lag2], lim );
end


% Integrate ...

Li = min([options.MaxLag floor(options.MaxLagFactor*LC)]);

% ... until limit crossing
if strcmp(options.Method,'integrate') 
    
    rho = nan(lag2,1);
    for lag = 0:lag2-1
        rho(lag+1) =  mean( xp(lag+1:Lx) .* yp(1:Lx-lag), 'omitnan' )/Fxy;
    end
    LS = trapz( [(0:lag2-1)';LC],[rho;0] );

% ... cumulatively until max lag and find maximum
elseif strcmp(options.Method,'cum-integrate')
    
    rho = nan(Li,1);
    for lag = 0:Li-1
        rho(lag+1) =  mean( xp(lag+1:Lx) .* yp(1:Lx-lag), 'omitnan' )/Fxy;
    end
    
    intR = cumtrapz( 0:Li-1, rho' );
    [~,LS] = max(intR);

% ... or take the limit crossing itself
else
    LS = LC;
end


% Multiply by DR to get length scale in physical units
LS = LS*options.dr;



% Plot

if options.Plot
    
    lw = 1.5;
    mks = 30;
    
    rv = (0:Li-1)*options.dr;
    
    [fig,~,co] = fig16x12('linlin',[1 1],'on','XLim',[0 rv(end)]);
    xlabel('$r$','Interpreter','latex')
    
    if ~strcmp(options.Method,'cum-integrate')  
        
        rho = nan(Li,1);
        for lag = 0:Li-1
            rho(lag+1) =  mean( xp(lag+1:Lx) .* yp(1:Lx-lag), 'omitnan' )/Fxy;
        end
      
        plot(rv,rho,'.','Color',co(1,:))
        plot([0,rv(end)],exp(-1)*[1 1],'LineStyle','--','LineWidth',lw,'Color','black')
        plot([0,rv(end)],[0 0],'LineStyle','--','LineWidth',lw,'Color','black')
        plot(LS,interp1(rv,rho,LS,'linear','extrap'),'.','MarkerSize',mks,'Color',co(2,:))
        
        ylabel('$\rho(r)$','Interpreter','latex')
        
    else
        
        plot(rv,intR*options.dr,'.','Color',co(1,:))
        plot(LS,max(intR)*options.dr,'.','MarkerSize',mks,'Color',co(2,:))
        
        ylabel('$\int_0^r\rho(r)$','Interpreter','latex')
        
    end
    
else
    fig = [];
end


end
