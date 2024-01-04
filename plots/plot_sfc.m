
function [fig,ax] = plot_sfc(S,U,vars,colors,varargin)

lw = 1;
ms = 4;

if nargin<4 || isempty(colors)
    colors = (1:100);
end


[fig,ax,co] = fig16x12('loglog',[1 0],varargin{:});

co = co(colors,:);

Nlvl = size(S,1);
Nvar = numel(vars);


for i_l = 1:Nlvl
    
    for i_v = 1:Nvar
        var = vars{i_v};
        c = co((i_l-1)*Nvar+i_v,:);
        
        y = S(i_l).(var);
        r = (1:length(y))*S(i_l).dr;
        plot(r(y>=0), y(y>=0),'o','Color',c,'MarkerSize',ms,'MarkerFaceColor',c,'HandleVisibility','off')
        plot(r(y<0), -y(y<0), 'o','Color',c,'MarkerSize',ms,'HandleVisibility','off')
        plot(nan,nan,'o','Color',c,'MarkerFaceColor',c)
        
        if ~isempty(U)
            u = U(i_l).(var);
            plot(r,abs(y)+u,'Color',c,'LineStyle','--','LineWidth',lw,'HandleVisibility','off')
            plot(r,abs(y)-u,'Color',c,'LineStyle','--','LineWidth',lw,'HandleVisibility','off')
        end
    end

end

xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')

end