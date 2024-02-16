
function [fig,ax] = plot_sfc(S,U,vars,colors,varargin)

lw = 1;
lwb = 1.5;
ms = 4;

if nargin<4 || isempty(colors)
    colors = (1:100);
end


[fig,ax,co] = fig16x12('loglog',[1 0],'on',varargin{:});

mk = repmat({'^','o','+','x','s'},1,5);

co = co(colors,:);
mk = mk(colors);

Nlvl = size(S,1);
Nvar = numel(vars);


for i_l = 1:Nlvl
    
    for i_v = 1:Nvar
        var = vars{i_v};
        if var(1)=='-'
            var = var(2:end);
            S(i_l).(var) = -S(i_l).(var);
        end
        c = co((i_l-1)*Nvar+i_v,:);
        
        y = S(i_l).(var);
        if isscalar(y)
            y = y*[1 1];
            r = ax.XLim;
            plot(r,y,'-','Color',c,'LineWidth',lwb)
        else
            r = (1:length(y))*S(i_l).dr;
            plot(r(y>=0), y(y>=0),'Color',c,'Marker',mk{i_v},'LineStyle','none',...
                'MarkerSize',ms,'MarkerFaceColor',c,'HandleVisibility','off')
            plot(r(y<0), -y(y<0),'Color',c,'Marker',mk{i_v},'LineStyle','none',...
                'MarkerSize',ms,'HandleVisibility','off')
            plot(nan,nan,'Color',c,'Marker',mk{i_v},'MarkerFaceColor',c,'LineStyle','none')
        end
        
        if ~isempty(U)
            u = U(i_l).(var);
            plot(r,abs(y)+u,'Color',c,'LineStyle','--','LineWidth',lw,'HandleVisibility','off')
            plot(r,abs(y)-u,'Color',c,'LineStyle','--','LineWidth',lw,'HandleVisibility','off')
        end
    end

end

xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')

end