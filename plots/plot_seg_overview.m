
function fig = plot_seg_overview(MOM,levels)

if nargin<2
    levels = unique(MOM.level);
end


fig = figure('Units','normalized','Position',[0 0 0.6 0.3]);
hold on, grid on
co = get(gca,'ColorOrder');

Nlvl = numel(levels);
for i_l = 1:Nlvl
    ind_l = (MOM.level==levels{i_l});
    plot(MOM{ind_l,"start"}, MOM{ind_l,"alt"},...
        'Marker','o','MarkerSize',10,'LineStyle','none',...
        'Color',co(i_l,:),'MarkerFaceColor',co(i_l,:))
    text(MOM{ind_l,"start"}, MOM{ind_l,"alt"}, MOM{ind_l,"name"})
end

legend(levels)
ylabel('Altitude [m]')

end