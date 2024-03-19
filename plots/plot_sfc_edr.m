
function [fig,ax] = plot_sfc_edr(S,U,vars,colors,npoints,varargin)

lwu = 1;
lwl = 1.5;
ms = 4;

if nargin<4 || isempty(colors)
    colors = (1:100);
end
if nargin<5 || isempty(npoints)
    npoints = 0;
end


[fig,ax,co] = fig16x12('loglog',[1 0],'on',varargin{:});

co = co(colors,:);

Nlvl = size(S,1);
Nvar = numel(vars);


for i_l = 1:Nlvl
    
    for i_v = 1:Nvar
        
        var = vars{i_v};
        c = co((i_l-1)*Nvar+i_v,:);
        
        fit_range = S(i_l).(strcat('fit_',var));
        slp       = S(i_l).(strcat('slp_',var));
        C         = S(i_l).(strcat('con_',var));
        edr       = S(i_l).(strcat('edr_',var));
        edr_free  = S(i_l).(strcat('edr_',var,'_free'));
        slp_free  = S(i_l).(strcat('slp_',var,'_free'));
        
        
        r = (1:length(S(i_l).(var)))*S(i_l).dr;
        ind_range = find( r>=fit_range(1) & r<=fit_range(2) );
        r_fit     = r(ind_range([1,end]));
                
        y = S(i_l).(var)./r.^slp;
        
        if npoints>0 && length(y)>1
            [isxlim,ind_xlim] = ismember('XLim',varargin(1:2:end));
            if isxlim
                rmin = varargin{2*ind_xlim}(1);
                rmax = varargin{2*ind_xlim}(2);
            else
                rmin = S(i_l).dr;
                rmax = length(y)*S(i_l).dr;
            end
            ind_r = unique(interp1(r,1:length(y),exp(linspace(log(rmin),log(rmax),npoints)),'nearest',length(y)));
            rs = r(ind_r);
            ys = y(ind_r);
        else
            rs = r;
            ys = y;
        end
        
        plot(rs,abs(ys),'o','MarkerSize',ms,'Color',c,'MarkerFaceColor',c)
        
        if ~isempty(U)
            u = U(i_l).(var)./r.^slp;
            plot(r,abs(y)+u,'Color',c,'LineStyle','--','LineWidth',lwu,'HandleVisibility','off')
            plot(r,abs(y)-u,'Color',c,'LineStyle','--','LineWidth',lwu,'HandleVisibility','off')
        end
        
        plot(r_fit,C*(r_fit*edr_free).^slp_free./r_fit.^slp, '--',...
            'LineWidth',lwl,'Color',c,'HandleVisibility','off')
        plot(r_fit,C*([1 1]*edr).^slp,'-',...
            'LineWidth',lwl,'Color',c,'HandleVisibility','off')

    end
    
end
    
xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')

end