
%% Settings

plot_vars = {'uu','vv','ww','wB','W','S3r'};


%% Load

addpath(genpath(myprojectpath))
load('S_eureca.mat')

plotpath = [myprojectpath,filesep,'figures',filesep,'errors'];
if ~isfolder(plotpath), mkdir(plotpath), end

Nseg = size(S,1);
Nvar = numel(plot_vars);
Nlvl = numel(levels);



%% Integrate / compensate

for i_s = 1:Nseg
    S(i_s).S3 = 3*( S(i_s).uuu + S(i_s).vvu + S(i_s).wwu );
    U(i_s).S3 = 3*sqrt( U(i_s).uuu.^2 + U(i_s).vvu.^2 + U(i_s).wwu.^2 );
end

for i_s = 1:Nseg
    r = (1:length(S(i_s).uuu))*S(i_s).dr;
    
    if MOM.RH_nan_fraction(i_s)==0
        S(i_s).W  = 6*cumtrapz(r,S(i_s).wB .*r.^2) ./ r.^3;
        S(i_s).Wt = 6*cumtrapz(r,S(i_s).wBt.*r.^2) ./ r.^3;
        S(i_s).Wq = 6*cumtrapz(r,S(i_s).wBq.*r.^2) ./ r.^3;
        U(i_s).W  = 6*cumtrapz(r,U(i_s).wB .*r.^2) ./ r.^3;
        U(i_s).Wt = 6*cumtrapz(r,U(i_s).wBt.*r.^2) ./ r.^3;
        U(i_s).Wq = 6*cumtrapz(r,U(i_s).wBq.*r.^2) ./ r.^3;
    end
       
    S(i_s).uuu3r = 3*S(i_s).uuu./r;
    S(i_s).vvu3r = 3*S(i_s).vvu./r;
    S(i_s).wwu3r = 3*S(i_s).wwu./r;
    U(i_s).uuu3r = 3*U(i_s).uuu./r;
    U(i_s).vvu3r = 3*U(i_s).vvu./r;
    U(i_s).wwu3r = 3*U(i_s).wwu./r;
    
    S(i_s).S3r = S(i_s).S3./r;
    U(i_s).S3r = U(i_s).S3./r;
end



%% PLOTS

mks = 2;
lw = 1;
Npoints = 40;

[f,~,co] = fig16x12;
close(f)
co = repmat(co,100,1);



%% Plot structure functions

for i_v = 1:Nvar
    var = plot_vars{i_v};
    
    for i_l = 1:Nlvl
        lvl = levels{i_l};
        
        fig = figure('Units','normalized','Position',[0 0 0.4 0.5]);
        grid on, hold on
        i_c = 1; legstr = [];
        
        for i_s = 1:Nseg
            if S(i_s).level==lvl
                c = co(i_c,:);
                plot(nan,nan,'o','Color',c,'MarkerFaceColor',c,'MarkerSize',mks)
                r = (1:length(S(i_s).uuu))*S(i_s).dr;
                y = S(i_s).(var);
                if ~isempty(y)
                    ind_r = unique(interp1(r,1:length(y),exp(linspace(log(r(1)),log(r(end)),Npoints)),'nearest',length(y)));
                    rs = r(ind_r); ys = y(ind_r);
                    plot(rs(ys>=0), ys(ys>=0),'o','Color',c,'LineStyle','none',...
                        'MarkerSize',mks,'MarkerFaceColor',c,'HandleVisibility','off')
                    plot(rs(ys<0),  -ys(ys<0),'o','Color',c,'LineStyle','none',...
                        'MarkerSize',mks,'HandleVisibility','off')
                end
                i_c = i_c + 1;
                legstr = vertcat(legstr,S(i_s).flight+" "+S(i_s).name);
            end
        end
        
        ind_l = find([avS(:).level]==lvl);
        i_c = 8; c = co(i_c,:);
        plot(nan,nan,'^','Color',c,'MarkerFaceColor',c,'MarkerSize',mks+3)
        r = (1:length(avS(ind_l).uuu))*avS(ind_l).dr;
        y = avS(ind_l).(var);
        ind_r = unique(interp1(r,1:length(y),exp(linspace(log(r(1)),log(r(end)),Npoints)),'nearest',length(y)));
        rs = r(ind_r); ys = y(ind_r);
        plot(rs(ys>=0), ys(ys>=0),'^','Color',c,'LineStyle','none',...
            'MarkerSize',mks+3,'MarkerFaceColor',c,'HandleVisibility','off')
        plot(rs(ys<0),  -ys(ys<0),'^','Color',c,'LineStyle','none',...
            'MarkerSize',mks+3,'HandleVisibility','off')
        legstr = vertcat(legstr,"AV");
        
        set(gca,'XScale','log','YScale','log')
        xlabel('r [m]')
        ylabel(['S ',var])
        legend(legstr,'Location','eastoutside')
        title(lvl)
        print(fig,[plotpath,filesep,'S_',var,'_',lvl],'-dpng','-r300')
    end
    
end



%% Integral length scale of increments

for i_v = 1:Nvar
    var = plot_vars{i_v};
    
    if isfield(L,var)
    
        for i_l = 1:Nlvl
            lvl = levels{i_l};

            fig = figure('Units','normalized','Position',[0 0 0.4 0.5]);
            grid on, hold on
            i_c = 1; legstr = [];

            for i_s = 1:Nseg
                if S(i_s).level==lvl
                    c = co(i_c,:);
                    r = (1:length(S(i_s).uuu))*S(i_s).dr;
                    if ~isempty(S(i_s).(var))
                        plot(r,L(i_s).(var),'Color',c,'LineWidth',lw);
                    else
                        plot(nan,nan,'.','Color',c);
                    end
                    i_c = i_c + 1;
                    legstr = vertcat(legstr,S(i_s).flight+" "+S(i_s).name);
                end
            end
            
            xlabel('r [m]')
            ylabel(['L ',var,' [m]'])
            legend(legstr,'Location','eastoutside')
            title(lvl)
            print(fig,[plotpath,filesep,'L_',var,'_',lvl],'-dpng','-r300')
        end
        
    end
    
end



%% Number of independent samples

for i_v = 1:Nvar
    var = plot_vars{i_v};
    
    if isfield(N,var)
    
        for i_l = 1:Nlvl
            lvl = levels{i_l};

            fig = figure('Units','normalized','Position',[0 0 0.4 0.5]);
            grid on, hold on
            i_c = 1; legstr = [];

            for i_s = 1:Nseg
                if S(i_s).level==lvl
                    c = co(i_c,:);
                    r = (1:length(S(i_s).uuu))*S(i_s).dr;
                    if ~isempty(S(i_s).(var))
                        plot(r,N(i_s).(var),'Color',c,'LineWidth',lw);
                    else
                        plot(nan,nan,'.','Color',c);
                    end
                    i_c = i_c + 1;
                    legstr = vertcat(legstr,S(i_s).flight+" "+S(i_s).name);
                end
            end
            
            ind_l = find([avS(:).level]==lvl);
            i_c = 8; c = co(i_c,:); 
            r = (1:length(avS(ind_l).uuu))*avS(ind_l).dr;
            plot(r,avN(ind_l).(var),'Color',c,'LineWidth',lw+1);
            legstr = vertcat(legstr,"AV");
            
            set(gca,'XScale','log','YScale','log')
            xlabel('r [m]')
            ylabel(['N ',var])
            legend(legstr,'Location','eastoutside')
            title(lvl)
            print(fig,[plotpath,filesep,'N_',var,'_',lvl],'-dpng','-r300')
        end
        
    end
    
end



%% Uncertainty

for i_v = 1:Nvar
    var = plot_vars{i_v};
    
    for i_l = 1:Nlvl
        lvl = levels{i_l};

        fig = figure('Units','normalized','Position',[0 0 0.4 0.5]);
        grid on, hold on
        i_c = 1; legstr = [];

        for i_s = 1:Nseg
            if S(i_s).level==lvl
                c = co(i_c,:);
                r = (1:length(S(i_s).uuu))*S(i_s).dr;
                if ~isempty(S(i_s).(var))
                    plot(r,U(i_s).(var),'Color',c,'LineWidth',lw);
                else
                    plot(nan,nan,'.','Color',c);
                end
                i_c = i_c + 1;
                legstr = vertcat(legstr,S(i_s).flight+" "+S(i_s).name);
            end
        end

        ind_l = find([avS(:).level]==lvl);
        i_c = 8; c = co(i_c,:); 
        r = (1:length(avS(ind_l).uuu))*avS(ind_l).dr;
        plot(r,avU(ind_l).(var),'-','Color',c,'LineWidth',lw+1);
        legstr = vertcat(legstr,"AV");

        set(gca,'XScale','log','YScale','log')
        xlabel('r [m]')
        ylabel(['U ',var])
        legend(legstr,'Location','eastoutside')
        title(lvl)
        print(fig,[plotpath,filesep,'U_',var,'_',lvl],'-dpng','-r300')
    end

end