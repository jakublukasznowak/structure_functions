
% MYPROJECTPATH is the path where you downloaded the codes

% Prepare paths

addpath(genpath(myprojectpath))

plotpath = [myprojectpath,filesep,'figures_Baltic'];
if ~isfolder(plotpath)
    mkdir(plotpath)
end



%% Load dataset

levels = {'level01'};

fsamp = 100;

load('BalticData1.mat')
TURB = BalticData';
load('BalticData2.mat')
TURB = cat(1,TURB,BalticData');
Level.Number = 1:numel(TURB);
clear BalticData


vars = fieldnames(TURB);

Nvar = numel(vars);
Nseg = numel(TURB);
Nlvl = numel(Level);


%% Label segments according to level

for i_l = 1:Nlvl
    for ii_s = 1:length(Level(i_l).Number)
        TURB(Level(i_l).Number(ii_s)).level = string(levels{i_l});
    end
end


%% Calculate segment averages

MOM = table([TURB(:).level]','VariableNames',{'level'});
for i_s = 1:Nseg
    for i_v = 1:Nvar
        var = vars{i_v};
        MOM.(var)(i_s) = mean(TURB(i_s).(var));
    end
    MOM.length(i_s) = length(TURB(i_s).U)*MOM.TAS(i_s)/fsamp;
end
MOM.dr = MOM.TAS/fsamp;
MOM.alt = MOM.Alt;


%% Calculate buoyancy

Mv = 18.015;     % [g/mol]
Md = 28.9647;
R  = 8.31432;     % [J/mol/K]
Rd = R/Md*1000;  % [J/kg/K]
Rv = R/Mv*1000;
ep = Rd/Rv;
cpd= 1003.6;    % [J/kg/K]
T0 = 273.15;     % [K]
g = 9.81;        % [m/s2]

for i_s = 1:Nseg  
    th = TURB(i_s).Theta+T0; % K
    q  = TURB(i_s).q; % kg/kg
        
    thv = th.*(1+(1/ep-1)*q); % K
    
    TURB(i_s).TH_V = thv; 

    TURB(i_s).B  = detrend( g/mean(thv)*(thv-mean(thv)) ); % buoyancy
    TURB(i_s).Bt = detrend( g/mean(th) *(th -mean(th)) );  % pseudobuoyancy t-part
    TURB(i_s).Bq = detrend( g/mean(th) *(1/ep-1)*(th.*q-mean(th.*q)) ); % pseudobuoyancy tq-part
end


% Detrend

for i_s = 1:Nseg
    TURB(i_s).U = detrend( TURB(i_s).U ) /3.6;
    TURB(i_s).V = detrend( TURB(i_s).V ) /3.6;
    TURB(i_s).W = detrend( TURB(i_s).W ) /3.6;
    TURB(i_s).p = detrend( TURB(i_s).p );
end



%% Calculate SFC with uncertainties for individual segments

% List of structure functions
sfc_vars = {'uuu',{'U','U','U'};
            'vvu',{'V','V','U'};
            'wwu',{'W','W','U'};
            'wB', {'W','B'};
            'wBt',{'W','Bt'};
            'wBq',{'W','Bq'};
            'uu', {'U','U'};
            'vv', {'V','V'};
            'ww', {'W','W'};
            };

% Max r separation
r_max = 1000;


S = table2struct(MOM(:,{'level','alt','dr','length'}));
V = S; L = S; N = S; U = S;

Nseg = size(S,1);
Nvar = size(sfc_vars,1);

us = etd(clock,0,Nseg*Nvar,30);

for i_v = 1:Nvar
    
    var = sfc_vars{i_v,1};
    varlist = sfc_vars{i_v,2};
    
    midFlag = startsWith(varlist,'X_');
    difFlag = ~startsWith(varlist,'X_');
    varlist(midFlag) = extractAfter(varlist(midFlag),'X_');
    
    fprintf('S_%s = < ',var)
    fprintf('%s_X ',varlist{midFlag})
    fprintf('d %s ',varlist{difFlag})
    fprintf('>\n')

    for i_s = 1:Nseg
        
        r_maxlag = ceil(r_max/MOM.dr(i_s))+1;
        
        Lt = length(TURB(i_s).U);
    
        A = cell2mat(cellfun(@(v) TURB(i_s).(v),varlist,'UniformOutput',false));
        
        if ~any(all(isnan(A),1))
            
            S(i_s).(var) = nan(1,r_maxlag);
            V(i_s).(var) = nan(1,r_maxlag);
            L(i_s).(var) = nan(1,r_maxlag);

            for r_lag = 1:r_maxlag
                
                difI = prod( A(r_lag+1:Lt,difFlag)-A(1:Lt-r_lag,difFlag) ,2);
                
                if sum(midFlag)>0
                    midI = prod( 0.5*(A(r_lag+1:Lt,midFlag)+A(1:Lt-r_lag,midFlag)) ,2);
                else
                    midI = 1;
                end
                
                I = difI.*midI; % increments
                
                S(i_s).(var)(r_lag) = mean(I,'omitnan');        % structure function
%                 V(i_s).(var)(r_lag) = mean(I.^2,'omitnan');     % variance of increments
%                 L(i_s).(var)(r_lag) = int_ls_short(I)*S(i_s).dr;% integral length scale for increments
            end
            
            % number of independent samples
            N(i_s).(var) = S(i_s).length./L(i_s).(var); 
            
            % uncertainty
            U(i_s).(var) = sqrt( (V(i_s).(var) - S(i_s).(var).^2) ./ N(i_s).(var) ); 
              
        end
        
        us = etd(us,Nseg*(i_v-1)+i_s);

    end
    
end


% Interpolation to common r-grid

dr = 0.4;
r = max(MOM.dr):dr:r_max;
r_maxlag = length(r);


for i_s = 1:Nseg
            
    for i_v = 1:Nvar
    
        var = sfc_vars{i_v,1};
        
        Lt = length(S(i_s).(var));

        r_old = (1:Lt)*S(i_s).dr;

        S(i_s).(var) = interp1( r_old, S(i_s).(var), r );
        V(i_s).(var) = interp1( r_old, V(i_s).(var), r );
        L(i_s).(var) = interp1( r_old, L(i_s).(var), r );
        
        N(i_s).(var) = S(i_s).length./L(i_s).(var); 
        U(i_s).(var) = sqrt( (V(i_s).(var) - S(i_s).(var).^2) ./ N(i_s).(var) ); 

    end
    
    S(i_s).dr = dr;
    
end



%% Average SFC at characteristic levels

Nseg = size(S,1);
Nlvl = numel(levels);
Nvar = size(sfc_vars,1);


avS = struct([]); avN = struct([]); avU = struct([]);

for i_l = 1:Nlvl
    ind_l = find([S(:).level]'==levels{i_l});
    
    avS(i_l,1).level  = string(levels{i_l});
    avS(i_l).number = numel(ind_l);
    avS(i_l).alt    = mean([S(ind_l).alt]);
    avS(i_l).dr     = mean([S(ind_l).dr]);
    avS(i_l).length = mean([S(ind_l).length]);
    
    for i_v = 1:Nvar
        var = sfc_vars{i_v,1};
        
        avS(i_l,1).(var) = mean( cat(1,S(ind_l).(var)), 1);
        avV              = mean( cat(1,V(ind_l).(var)), 1);
        avN(i_l,1).(var) = sum( repmat([S(ind_l).length]',1,r_maxlag) ...
            ./ cat(1,L(ind_l).(var)) );
        avU(i_l,1).(var) = sqrt( (avV - avS(i_l).(var).^2) ./ avN(i_l).(var) );
    end
    
    % s3L
    if all(ismember({'uuu','vvu','wwu'},fieldnames(avS)))
        avS(i_l).s3l = 3*( avS(i_l).uuu + avS(i_l).vvu + avS(i_l).wwu );
        avU(i_l).s3l = 3*sqrt( avU(i_l).uuu.^2 + avU(i_l).vvu.^2 + avU(i_l).wwu.^2 );
    end
    
    % wBt, wBq
    if all(ismember({'wBt','wBq'},fieldnames(avS)))
        avS(i_l).wQ  = avS(i_l).wBt + avS(i_l).wBq;
        avU(i_l).wQ  = sqrt(avU(i_l).wBt.^2 + avU(i_l).wBq.^2);
    end
    
    % transport
    if all(ismember({'uuW','vvW','wwW'},fieldnames(avS)))
        avS(i_l).s2W  = avS(i_l).uuW + avS(i_l).vvW + avS(i_l).wwW;
        avU(i_l).s2W  = sqrt(avU(i_l).uuW.^2 + avU(i_l).vvW.^2 + avU(i_l).wwW.^2);
    end
end



%% Compensate / integrate

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    
    avS(i_l).W  = 6*cumtrapz(r,avS(i_l).wB .*r.^2) ./ r.^3;
    avS(i_l).Wt = 6*cumtrapz(r,avS(i_l).wBt.*r.^2) ./ r.^3;
    avS(i_l).Wq = 6*cumtrapz(r,avS(i_l).wBq.*r.^2) ./ r.^3;
    avU(i_l).W  = 6*cumtrapz(r,avU(i_l).wB .*r.^2) ./ r.^3;
    avU(i_l).Wt = 6*cumtrapz(r,avU(i_l).wBt.*r.^2) ./ r.^3;
    avU(i_l).Wq = 6*cumtrapz(r,avU(i_l).wBq.*r.^2) ./ r.^3;
    
    avS(i_l).s3lr = avS(i_l).s3l./r;
    avU(i_l).s3lr = avU(i_l).s3l./r;
    
    avS(i_l).Ws3lr = avS(i_l).W - avS(i_l).s3lr;
    avU(i_l).Ws3lr = avU(i_l).W + avU(i_l).s3lr;
      
    avS(i_l).Wrs3l = avS(i_l).W.*r - avS(i_l).s3l;
    avU(i_l).Wrs3l = avU(i_l).W.*r + avU(i_l).s3l;
    
end



%% Calculate dissipation rate

fit_range = [10 30];

edr_vars = {'uu','vv','ww','Wrs3l'};
slps = [2/3 2/3 2/3 1];
cons = [2.0 2.6 2.6 4.0];

Nvar = numel(edr_vars);
for i_v = 1:Nvar
    var = edr_vars{i_v};
    [avS,avU] = calc_edr(avS,avU,edr_vars{i_v},'Slope',slps(i_v),'Constant',cons(i_v),...
        'FitRange',fit_range,'Method','direct','FitPoints',6);
end


% Compensate for plotting
Nlvl = size(avS,1);
for i_l = 1:Nlvl
    for i_v = 1:Nvar
        var = edr_vars{i_v};
        avS(i_l).([var,'_c']) = avS(i_l).(var) ./ r.^slps(i_v) / cons(i_v);
        avU(i_l).([var,'_c']) = avU(i_l).(var) ./ r.^slps(i_v) / cons(i_v);
    end
    avS(i_l).edr_s2_4 = avS(i_l).edr_s2*4;
    avU(i_l).edr_s2_4 = avU(i_l).edr_s2*4;
end



%% Save/load

% save('S_baltic_1km.mat','S','U','N','L','V','avS','avU','avN',...
%     'r_max','dr','r_maxlag','r','MOM','levels','plotpath')

load('S_baltic_1km.mat')
addpath(genpath(myprojectpath))
plotpath = [myprojectpath,filesep,'figures_Baltic'];
if ~isfolder(plotpath)
    mkdir(plotpath)
end


%% TABLES

% Segment overview

print_table(MOM,{'length_km','alt'},1,0)


% Dissipation rates

print_vars = {'uu','vv','ww','s2','Wrs3l'};

Nlvl = size(avS,1);
Nvar = numel(print_vars);

fprintf(' %20s',''), fprintf(' & %14s',print_vars{:}), fprintf(' \\\\ \n')
for i_l = 1:Nlvl
    fprintf(' %20s',avS(i_l).level)
    for i_v = 1:Nvar
        var = print_vars{i_v};
        fprintf(' & %6.*f (%.*f)',3,1e4*avS(i_l).(['edr_',var]),3,1e4*avU(i_l).(['edr_',var]) )
    end
    fprintf(' \\\\ \n')
end



%% PLOTS


%% Structure functions

xlim = [dr r_max];
ylim = [1e-8 1e-1];

Npoints = 40;

Nlvl = size(avS,1);


for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','Wt','Wq'},[7 5 1],Npoints,'XLim',xlim,'YLim',ylim);
    legend({'$W$','$W_\theta$','$W_q$'},'Interpreter','latex','Location','northwest')
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'W_',levels{i_l}],'-dpng','-r300')
end


for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'s3lr'},1,Npoints,'XLim',xlim,'YLim',ylim);
    legend({'$S_3 r^{-1}$'},'Interpreter','latex','Location','southeast')
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s3l_',levels{i_l}],'-dpng','-r300')
end


for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'Ws3lr','edr_s2_4'},[7 3],Npoints,'XLim',xlim,'YLim',ylim);
    legend({'$W-S_3 r^{-1}$','$4\epsilon_2$'},'Interpreter','latex','Location','southeast')
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'Ws3l_',levels{i_l}],'-dpng','-r300')
end


for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','-s3lr','edr_s2_4'},[7 1 3],Npoints,'XLim',xlim,'YLim',ylim);
    legend({'$W$','$-S_3 r^{-1}$','$4\epsilon_2$'},'Interpreter','latex','Location','southeast')
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'balance_',levels{i_l}],'-dpng','-r300')
end


for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'uu_c','vv_c','ww_c'},[7 5 1],Npoints,'XLim',xlim,'YLim',[1e-3 4e-2]);
    legend({'uu','vv','ww'},'Interpreter','latex','Location','southeast')
    ylabel('$S_2r^{-2/3}C^{-1}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s2_',levels{i_l}],'-dpng','-r300')
end




%% Fitting dissipation

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    
    fig = plot_sfc_edr(avS(i_l),[],{'uu','vv','ww'},[7 5 1],Npoints,'XLim',xlim);
    legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    ylabel('$S_2r^{-2/3}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr2_',levels{i_l}],'-dpng','-r300')
    
%     fig = plot_sfc_edr(avS(i_l),avU(i_l),{'Wrs3l'},7,Npoints,'XLim',xlim);
%     ylabel('$W-S_3^L r^{-1}$','Interpreter','latex')
%     title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
%     print(fig,[plotpath,filesep,'edr3_',levels{i_l}],'-dpng','-r300')
    
end

