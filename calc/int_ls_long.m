
function [LS,fig] = int_ls_long (rho,options)

arguments
    rho (:,1) {mustBeReal, mustBeNonempty}
    options.Method (1,1) string {mustBeMember(options.Method,...
        {'e-decay','zero','integrate','cum-integrate'})} = 'e-decay'
    options.MaxLag (1,1) {mustBeInteger, mustBeFinite} = length(rho)-1
    options.dr (1,1) {mustBePositive} = 1
    options.Plot (1,1) logical = false
end


lagvec = (0:options.MaxLag)';


% Select critical limit value

if strcmp(options.Method,'e-decay')
    lim = exp(-1);
else
    lim = 0;
end


% Find the crossing

ind0 = find(rho<lim,1,'first');

if isempty(ind0)
    LC = options.MaxLag;
    
    if ismember(options.Method,{'e-decay','zero'})
        warning('INT_LS_LONG:NoCrossing','\nNo %.1f crossing. Take max lag %.1f.',lim,LC)
    elseif strcmp(options.Method,'integrate')
        warning('INT_LS_LONG:NoCrossing','\nNo %.1f crossing. Integrating to max lag %.1f.',lim,LC)
    end
    
else
    LC = interp1( rho(ind0-1:ind0), lagvec(ind0-1:ind0), lim );
end


% Integrate ...

% ... until limit crossing
if strcmp(options.Method,'integrate') 
    
    LS = trapz( [lagvec(1:ind0-1);LC],[rho(1:ind0-1);lim] );

% ... cumulatively until max lag and find maximum
elseif strcmp(options.Method,'cum-integrate')
    
    intR = cumtrapz( lagvec, rho );
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
    
    rv = lagvec*options.dr;
    
    [fig,~,co] = fig16x12('linlin',[1 1],'on','XLim',[0 rv(end)]);
    xlabel('$r$','Interpreter','latex')
    
    if ~strcmp(options.Method,'cum-integrate')  

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