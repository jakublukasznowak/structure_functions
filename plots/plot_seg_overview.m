
function fig = plot_seg_overview(MOM,levels,ifdirs,ifnames)

if nargin<4 || isempty(ifnames)
    ifnames = false;
end
if nargin<3 || isempty(ifdirs)
    ifdirs = true;
end
if nargin<2
    levels = unique(MOM.level);
end


mk = 'o';
mks = 10;

dirs   = {'along','cross'};


[fig,ax,co] = fig16x12;
fig.PaperSize = [36 10]; fig.PaperPosition = [0 0 36 10];
ax.Position = [0.07 0.14 1-1.5*0.07 1-1.5*0.14];


Nlvl = numel(levels);
dt = datetime(clock,'TimeZone','UTC');
if Nlvl>1
    for i_l = 1:Nlvl
        c = co(i_l,:);
        plot(dt,nan,'Color',c,'MarkerFaceColor',c,...
            'Marker',mk,'MarkerSize',mks,'LineStyle','none')
    end
end
if ifdirs
    plot(dt,nan,'Color','black','MarkerFaceColor','black',...
        'Marker',mk,'MarkerSize',mks,'LineStyle','none')
    plot(dt,nan,'Color','black','MarkerFaceColor','none',...
        'Marker',mk,'MarkerSize',mks,'LineStyle','none')
end


for i_l = 1:Nlvl
    c = co(i_l,:);
    if ifdirs
        ind = MOM.level==levels{i_l} & MOM.dir2==dirs{1};
        plot(MOM.start(ind), MOM.alt(ind),'Marker',mk,'MarkerSize',mks,...
            'LineStyle','none','Color',c,'MarkerFaceColor',c)
        
        ind = MOM.level==levels{i_l} & MOM.dir2==dirs{2};
        plot(MOM.start(ind), MOM.alt(ind),'Marker',mk,'MarkerSize',mks,...
            'LineStyle','none','Color',c,'MarkerFaceColor','none')
    else
        ind = MOM.level==levels{i_l};
        plot(MOM.start(ind), MOM.alt(ind),'Marker',mk,'MarkerSize',mks,...
            'LineStyle','none','Color',c,'MarkerFaceColor',c)
    end
    
    ind = MOM.level==levels{i_l};
    if ifnames && ismember('name',MOM.Properties.VariableNames)
        text(MOM.start(ind), MOM.alt(ind), MOM.name(ind))
    end
end

legend(levels,'Location','best')
ylabel('Altitude [m]')

end