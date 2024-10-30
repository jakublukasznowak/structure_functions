

addpath(genpath(myprojectpath))
plotpath = [myprojectpath,filesep,'figures_bomex'];
if ~isfolder(plotpath), mkdir(plotpath), end


% Load profiles

p = load('bomex_profiles.mat');
Nt = size(p.thl,2);


% Constants

Rd = 287.04;  % [J/kg/K]
Rv = 461.5;
ep = Rd/Rv;
cpd= 1005;    % [J/kg/K]
Lv = 2.5e6;   % [J/kg]
T0 = 273.15;  % [K]
g  = 9.81;    % [m/s2]
p_ref = 1000e2; % [Pa]



%% Estimate characteristic levels

% Thresholds

stb_thresh = 0.1;     % [K/hPa]

% qv_thresh  = 0.35e-3; % [kg/kg]
% th_thresh  = 0.15/2;  % [K]
% thv_thresh = 0.2/2;   % [K]

iz_ml = (p.z>=100 & p.z<=300);
iz_cl = (p.z>=1000 & p.z<=1200);

ml_qv = mean(p.qv(iz_ml,:)); ml_th = mean(p.th(iz_ml,:)); ml_thv = mean(p.thv(iz_ml,:));
cl_qv = mean(p.qv(iz_cl,:)); cl_th = mean(p.th(iz_cl,:)); cl_thv = mean(p.thv(iz_cl,:));

% qv_thresh  = mean( 0.1*abs(ml_qv-cl_qv) );
% thv_thresh = mean( 0.1*abs(ml_thv-cl_thv) );
qv_thresh  = 0.1*abs(ml_qv(end)-cl_qv(end));
thv_thresh = 0.1*abs(ml_thv(end)-cl_thv(end));
th_thresh  = thv_thresh;


% CB: cloud base
    
% (1) Max cloud fraction
[cfmax,p.cb_cf_iz] = max(p.cf);
p.cb_cf = p.z(p.cb_cf_iz)';
p.cb_cf(cfmax==0) = nan;

% (2) Lifting condensation level for the conditions between 50 and 300 m
% Albright et al. 2022: sec. 3(a); Bony et al. 2022: Tab. 3
p.T_lcl = 1 ./ ( 1./(p.T-55) - (log(p.rh)/2840) ) + 55; % Bolton 1980
p.z_lcl = repmat(p.z,1,Nt) + cpd/g*(p.T-p.T_lcl);
iz = find(p.z>=50 & p.z<=300);
p.cb_lcl = mean(p.z_lcl(iz,:));


% ML: mixed layer top

iz100 = find(p.z>=100,1,'first');

% (1) Gradient method for specific humidity
% Albright et al. 2022: sec. 3(a) and A(a)
qvm = cumtrapz(p.z(iz100:end),p.rho(iz100:end,:).*p.qv(iz100:end,:)) ....
    ./ cumtrapz(p.z(iz100:end),p.rho(iz100:end,:));
qvm = [p.qv(1:iz100,:);qvm(2:end,:)];
[~,p.ml_qv_iz] = max(abs(p.qv-qvm)>qv_thresh);
p.ml_qv = p.z(p.ml_qv_iz)';

% (2) Gradient method for potential temperature
% Albright et al. 2022: sec. 3(a) and A(a)
thm = cumtrapz(p.z(iz100:end),p.rho(iz100:end,:).*p.th(iz100:end,:)) ....
    ./ cumtrapz(p.z(iz100:end),p.rho(iz100:end,:));
thm = [p.th(1:iz100,:);thm(2:end,:)];
[~,p.ml_th_iz] = max(abs(p.th-thm)>th_thresh);
p.ml_th = p.z(p.ml_th_iz)';

% (3) Max relative humidity
% c.f. Bony et al. 2022: fig. 3; Albright et al. 2022: sec. A(c)
[~,p.ml_rh_iz] = max([zeros(iz100-1,Nt);p.rh(iz100:end,:)]);
p.ml_rh = p.z(p.ml_rh_iz)';


% INV: trade-wind inversion

% (1) Strong enough static stability
% c.f. Albright et al. 2022: sec. A(b)
p.stb = -diff(p.thl)./diff(p.p/100);
[~,p.inv_stb_iz] = max(p.stb>stb_thresh);
p.inv_stb = 0.5*(p.z(p.inv_stb_iz)+p.z(p.inv_stb_iz+1))';

% (2) Min moist static energy between 1400 and 4000 m
% Bony et al. 2022: Tab. 3
p.s = cpd*p.T+g*repmat(p.z,1,Nt)+Lv*p.qv;
iz1300 = find(p.z>=1300,1,'first');
iz4000 = find(p.z<=4000,1,'last');
[~,p.inv_s_iz] = min([inf*ones(iz1300-1,Nt);p.s(iz1300:iz4000,:)]);
p.inv_s = p.z(p.inv_s_iz)';

% (3) Max gradient of virtual potential temperature
% Bony et al. 2022: fig. 3
[~,p.inv_thv_iz] = max(diff(p.thv)./diff(p.z));
p.inv_thv = 0.5*(p.z(p.inv_thv_iz)+p.z(p.inv_thv_iz+1))';

% (4) Max gradient of relative humidity
% Bony et al. 2022: fig. 3
[~,p.inv_rh_iz] = max(-diff(p.rh)./diff(p.z));
p.inv_rh = 0.5*(p.z(p.inv_rh_iz)+p.z(p.inv_rh_iz+1))';


% SC: subcloud layer top

% (1) Gradient method for virtual potential temperature
% Albright et al. 2022: sec. 3(a) and A(a)
thvm = cumtrapz(p.z(iz100:end),p.rho(iz100:end,:).*p.thv(iz100:end,:)) ....
    ./ cumtrapz(p.z(iz100:end),p.rho(iz100:end,:));
thvm = [p.thv(1:iz100,:);thvm(2:end,:)];
[~,p.sc_thv_iz] = max(abs(p.thv-thvm)>thv_thresh);
p.sc_thv = p.z(p.sc_thv_iz)';

% (2) Level of neutral buoyancy
% Albright et al. 2022: sec. 3(a) and A(b)
iz50 = find(p.z>=50,1,'first');
thv_sfc = mean(p.thv(1:iz50,:));
p.sc_nb = nan(1,Nt);
for it = 1:Nt
    iz = find( p.z>p.ml_rh(it) & p.z<=p.inv_stb(it) );
    pp = polyfit(p.z(iz),p.thv(iz,it),1);
    p.sc_nb(it) = (thv_sfc(it)-pp(2))/pp(1);
end


save('bomex_profiles.mat','-struct','p')



%% Plot profiles

for it = 1:Nt

    fig = figure('Color','white','PaperUnits','centimeters',...
        'PaperSize',[1.5*16 1.5*9],'PaperPosition',[0 0 1.5*16 1.5*9]);
    t = tiledlayout(fig,1,5);

    set(groot,'defaultLineLineWidth',2);

    ax1 = nexttile; hold on, grid on, box on
    co = get(ax1,'ColorOrder');
    plot(ax1,p.thl(:,it)-T0,p.z,'Color',co(2,:))
    plot(ax1,p.thv(:,it)-T0,p.z,'Color',co(3,:))
    ylim([0 2200]), xlim([25.5 36])
    xlabel('$\theta_l,\,\theta_v\,[\mathrm{^\circ C}]$','Interpreter','latex')

    ax2 = nexttile; hold on, grid on, box on
    plot(ax2,1e3*p.qv(:,it),p.z,'Color',co(4,:))
    ylim(ax1.YLim), xlim([4 18])
    yticklabels(ax2,{})
    xlabel('$q_v\,[\mathrm{g/kg}]$','Interpreter','latex')

    ax3 = nexttile; hold on, grid on, box on
    plot(ax3,1e6*p.ql(:,it),p.z,'Color',co(1,:))
    ylim(ax1.YLim), xlim([0 20])
    yticklabels(ax3,{})
    xlabel('$q_l\,[\mathrm{g/t}]$','Interpreter','latex')

    ax4 = nexttile; hold on, grid on, box on
    plot(ax4,100*p.rh(:,it),p.z,'Color',co(4,:))
    ylim(ax1.YLim), xlim([29 100])
    yticklabels(ax4,{})
    xlabel('$RH\,[\mathrm{g/t}]$','Interpreter','latex')

    ax5 = nexttile; hold on, grid on, box on
    plot(ax5,100*p.cf(:,it),p.z,'Color',co(1,:))
    ylim(ax1.YLim), xlim([0 10])
    yticklabels(ax5,{})
    xlabel('$CF\,[\mathrm{\%}]$','Interpreter','latex')
    
    xl1=ax1.XLim; xl2=ax2.XLim; xl3=ax3.XLim; xl4=ax4.XLim; xl5=ax5.XLim;
    
    plot(ax2,ax2.XLim,[1 1]*p.ml_qv(it),'-','Color',co(2,:))
    plot(ax1,ax1.XLim,[1 1]*p.ml_th(it),'--','Color',co(2,:))
    plot(ax4,ax4.XLim,[1 1]*p.ml_rh(it),':','Color',co(2,:))

    plot(ax1,ax1.XLim,[1 1]*p.sc_thv(it),'-','Color',co(3,:))
    plot(ax1,ax1.XLim,[1 1]*p.sc_nb(it),'--','Color',co(3,:))
    
    plot(ax5,ax5.XLim,[1 1]*p.cb_cf(it),'-','Color',co(1,:))
    plot(ax3,ax3.XLim,[1 1]*p.cb_lcl(it),'--','Color',co(1,:))

    plot(ax1,ax1.XLim,[1 1]*p.inv_stb(it),'-','Color',co(4,:))
    plot(ax2,ax2.XLim,[1 1]*p.inv_s(it),'--','Color',co(4,:))
    plot(ax1,ax1.XLim,[1 1]*p.inv_thv(it),':','Color',co(4,:))
    plot(ax4,ax4.XLim,[1 1]*p.inv_rh(it),'-.','Color',co(4,:))

    ax1.XLim=xl1; ax2.XLim=xl2; ax3.XLim=xl3; ax4.XLim=xl4; ax5.XLim=xl5;
    
    linkaxes([ax1,ax2,ax3,ax4,ax5],'y')

    title(t,sprintf('t = %.2f h',p.time(it)/3600),'Interpreter','latex')
    t.TileSpacing = 'none';
    t.Padding = 'none';

    print([plotpath,filesep,sprintf('prof_%03dmin',p.time(it)/60)],'-dpng','-r300')
end



%% Plot layers

set(groot,'defaultLineLineWidth',2);

fig = figure('Color','white','PaperUnits','centimeters',...
        'PaperSize',[1.5*16 1.5*9],'PaperPosition',[0 0 1.5*16 1.5*9]);
ax = axes('Position',[0.05 0.08 1-1.5*0.05 1-1.5*0.08]);
hold on, grid on, box on
co = get(gca,'ColorOrder');


time = p.time/3600;

plot(time,p.ml_qv,'-o','Color',co(2,:))
plot(time,p.ml_th,'--x','Color',co(2,:))
plot(time,p.ml_rh,':+','Color',co(2,:))

plot(time,p.sc_thv,'-o','Color',co(3,:))
plot(time,p.sc_nb,'--x','Color',co(3,:))

plot(time,p.cb_cf,'-o','Color',co(1,:))
plot(time,p.cb_lcl,'--x','Color',co(1,:))

% plot(time,p.inv_stb,'-o','Color',co(4,:))
% plot(time,p.inv_s,'--x','Color',co(4,:))
% plot(time,p.inv_thv,':+','Color',co(4,:))
% plot(time,p.inv_rh,'-.^','Color',co(4,:))

legend({'ML:$q_v$','ML:$\theta$','ML:RH', 'SC:$\theta_v$','SC:NB',...
    'CB:CF','CB:LCL', 'INV:stb','INV:s','INV:$\theta_v$','INV:RH'},...
    'Interpreter','latex','Location','best')
xlabel('Time [h]')

print([plotpath,filesep,'layers'],'-dpng','-r300')




%%

% set(groot,'defaultLineMarkerSize',4);
% 
% [~,~,co] = fig16x12;
% plot(h.thl,h.z,'o','Color',co(1,:))
% plot(h.thv,h.z,'o','Color',co(2,:))
% plot(h.th, h.z,'o','Color',co(3,:))
% plot(h.T,  h.z,'o','Color',co(4,:))
% plot(h.Tv, h.z,'o','Color',co(5,:))
% plot(f.thl,f.z,'o','Color',co(1,:),'MarkerFaceColor',co(1,:))
% plot(f.thv,f.z,'o','Color',co(2,:),'MarkerFaceColor',co(2,:))
% plot(f.th, f.z,'o','Color',co(3,:),'MarkerFaceColor',co(3,:))
% plot(f.T  ,f.z,'o','Color',co(4,:),'MarkerFaceColor',co(4,:))
% plot(f.Tv ,f.z,'o','Color',co(5,:),'MarkerFaceColor',co(5,:))
% legend({'thl','thv','th','T','Tv'})
% 
% [~,~,co] = fig16x12;
% plot(1e3*h.qt,h.z,'o','Color',co(4,:))
% plot(1e6*h.ql,h.z,'o','Color',co(1,:))
% plot(1e3*f.qt,f.z,'o','Color',co(4,:),'MarkerFaceColor',co(4,:))
% plot(1e6*f.ql,f.z,'o','Color',co(1,:),'MarkerFaceColor',co(1,:))
% legend({'qt','1e3*ql'})
% 
% [~,~,co] = fig16x12;
% plot(100*h.rh,h.z,'o','Color',co(4,:))
% plot(100*f.rh,f.z,'o','Color',co(4,:),'MarkerFaceColor',co(4,:))
% plot(100*h.cf,h.z,'o','Color',co(1,:))
% plot(100*f.cf,f.z,'o','Color',co(1,:),'MarkerFaceColor',co(1,:))