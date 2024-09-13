

%% Settings

% Compute uncertainties?
% (takes quite a long time)
if_uncertainties = true;

% Structure functions
r_max = 1000; % [m] Max r separation
dr = 0.4;     % [m] Interpolation spacing

% Dissipation rate
edr_fit_range = [10 30]; % [m]

% Sampling rate
fsamp = 100;



%% Prepare paths

addpath(genpath(myprojectpath))

datapath = mydatapath;

plotpath = [myprojectpath,filesep,'figures_baltic'];
if ~isfolder(plotpath), mkdir(plotpath), end



%% Load dataset

load('BalticData1.mat')
TURB = BalticData';
load('BalticData2.mat')
TURB = vertcat(TURB,BalticData');

Level.Number = 1:numel(TURB);


% Label segments according to level

levels = {'level01'};

Nlvl = numel(Level);

for i_l = 1:Nlvl
    for ii_s = 1:length(Level(i_l).Number)
        TURB(Level(i_l).Number(ii_s)).level = string(levels{i_l});
    end
end


clear BalticData Level


% Detrend

Nseg = numel(TURB);

for i_s = 1:Nseg
    TURB(i_s).U = detrend( TURB(i_s).U );
    TURB(i_s).V = detrend( TURB(i_s).V );
    TURB(i_s).W = detrend( TURB(i_s).W );
    TURB(i_s).p = detrend( TURB(i_s).p );
end



%% Compute auxilliary parameters

% Segment averages

MOM = table([TURB(:).level]','VariableNames',{'level'});

vars = fieldnames(TURB);

Nvar = numel(vars);
Nseg = numel(TURB);

for i_s = 1:Nseg
    for i_v = 1:Nvar
        var = vars{i_v};
        if numel(TURB(i_s).(var))>1
            MOM.(var)(i_s) = mean(TURB(i_s).(var));
        else
            MOM.(var)(i_s) = TURB(i_s).(var);
        end
    end
    MOM.length(i_s) = length(TURB(i_s).U)*MOM.TAS(i_s)/fsamp;
end

MOM.dr = MOM.TAS/fsamp;
MOM.alt = MOM.Alt;


% Buoyancy

R  = 8.31432;                   % [J/mol/K]
Mv = 18.015; Md = 28.9647;      % [g/mol]
Rd = R/Md*1000; Rv = R/Mv*1000; % [J/kg/K]
ep = Rd/Rv;
cpd= 1003.6; % [J/kg/K]
T0 = 273.15; % [K]
g  = 9.81;   % [m/s2]

for i_s = 1:Nseg  
    th = TURB(i_s).Theta+T0; % K
    q  = TURB(i_s).q; % kg/kg
        
    thv = th.*(1+(1/ep-1)*q); % K
    
    TURB(i_s).THV = thv; 

    TURB(i_s).B  = detrend( g/mean(thv)*(thv-mean(thv)) ); % buoyancy
    TURB(i_s).Bt = detrend( g/mean(thv) *(th -mean(th)) );  % pseudobuoyancy t-part
    TURB(i_s).Bq = detrend( g/mean(thv) *(1/ep-1)*(th.*q-mean(th.*q)) ); % pseudobuoyancy tq-part
end



%% Compute structure functions for individual segments

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
            'wXuu',{'U','U','X_W'};
            'wXvv',{'V','V','X_W'};
            'wXww',{'W','W','X_W'}};

        
% Prepare data structures

S = table2struct(MOM(:,{'level','alt','dr','length'}));
V = S; L = S; N = S; U = S;

% Compute

Nseg = size(S,1);
Nvar = size(sfc_vars,1);
us = etd(clock,0,Nseg*Nvar,30);

for i_v = 1:Nvar
    
    var = sfc_vars{i_v,1};
    varlist = sfc_vars{i_v,2};
    
    midFlag = startsWith(varlist,'X_');  % factors of type w_X
    absFlag = startsWith(varlist,'A_');  % factors of type |dw|
    difFlag = ~(midFlag | absFlag);      % factors of type dw
    varlist(midFlag) = extractAfter(varlist(midFlag),'X_');
    varlist(absFlag) = extractAfter(varlist(absFlag),'A_');
    
    fprintf('S_%s = < ',var)
    if sum(midFlag), fprintf('%s_X ',  varlist{midFlag}), end
    if sum(absFlag), fprintf('|d %s| ',varlist{absFlag}), end
    if sum(difFlag), fprintf('d %s ',  varlist{difFlag}), end
    fprintf('>\n')

    for i_s = 1:Nseg
        
        r_maxlag = ceil(r_max/MOM.dr(i_s))+1;
        
        A = cell2mat(cellfun(@(v) TURB(i_s).(v),varlist,'UniformOutput',false));

        Lt = length(TURB(i_s).U);
        
        if ~any(all(isnan(A),1))
            
            S(i_s).(var) = nan(1,r_maxlag);
            V(i_s).(var) = nan(1,r_maxlag);
            L(i_s).(var) = nan(1,r_maxlag);
            S(i_s).N     = nan(1,r_maxlag);

            for r_lag = 1:r_maxlag
                
                midI = 1; absI = 1; difI = 1;
                  
                % factors of type w_X
                if sum(midFlag)>0
                    midI = prod( 0.5*(A(r_lag+1:Lt,midFlag)+A(1:Lt-r_lag,midFlag)) ,2);
                end
                
                % factors of type |dw|
                if sum(absFlag)>0
                    absI = prod( abs( A(r_lag+1:Lt,absFlag)-A(1:Lt-r_lag,absFlag) ) ,2);
                end
                
                % factors of type dw
                if sum(difFlag)>0
                    difI = prod( A(r_lag+1:Lt,difFlag)-A(1:Lt-r_lag,difFlag) ,2);
                end
                
                % increments
                I = midI.*absI.*difI; 
                
                % structure function
                S(i_s).(var)(r_lag) = mean(I,'omitnan');        
                
                % parameters needed for uncertainty
                if if_uncertainties
                    V(i_s).(var)(r_lag) = mean(I.^2,'omitnan');      % variance of increments
                    L(i_s).(var)(r_lag) = int_ls_short(I)*S(i_s).dr; % integral length scale for increments
                    S(i_s).N    (r_lag) = sum(~isnan(I));            % number of valid samples
                end
                
            end
            
            % number of independent samples
            N(i_s).(var) = S(i_s).length./L(i_s).(var); 
            
            % uncertainty
            U(i_s).(var) = sqrt( (V(i_s).(var) - S(i_s).(var).^2) ./ N(i_s).(var) ); 
              
        end
        
        us = etd(us,Nseg*(i_v-1)+i_s);

    end
    
end


% Interpolate to common r-grid

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



%% Average at characteristic levels

avS = struct([]); avN = struct([]); avU = struct([]);

Nlvl = numel(levels);
Nvar = size(sfc_vars,1);

for i_l = 1:Nlvl
    ind_l = find([S(:).level]'==levels{i_l});
    
    avS(i_l,1).level = string(levels{i_l});
    avS(i_l).number  = numel(ind_l);
    avS(i_l).alt     = mean([S(ind_l).alt]);
    avS(i_l).dr      = mean([S(ind_l).dr]);
    avS(i_l).length  = mean([S(ind_l).length]);
    
    for i_v = 1:Nvar
        var = sfc_vars{i_v,1};
        ind_lv = ind_l;
        
        avS(i_l,1).(var) = mean( vertcat(S(ind_lv).(var)), 1);
        avV              = mean( vertcat(V(ind_lv).(var)), 1);
%         avN(i_l,1).(var) = sum( vertcat(N(ind_lv).(var)) ./ vertcat(L(ind_lv).(var)) );
        avN(i_l,1).(var) = sum( repmat([S(ind_lv).length]',1,r_maxlag) ./ vertcat(L(ind_lv).(var)) );
        avU(i_l,1).(var) = sqrt( (avV - avS(i_l).(var).^2) ./ avN(i_l).(var) );
    end
    
    % s3
    if all(ismember({'uuu','vvu','wwu'},fieldnames(avS)))
        avS(i_l).S3 = 3*( avS(i_l).uuu + avS(i_l).vvu + avS(i_l).wwu );
        avU(i_l).S3 = 3*sqrt( avU(i_l).uuu.^2 + avU(i_l).vvu.^2 + avU(i_l).wwu.^2 );
    end
    
    % Transport integrand
    if all(ismember({'wXuu','wXvv','wXww'},fieldnames(avS)))
        avS(i_l).wXS2  = avS(i_l).wXuu + avS(i_l).wXvv + avS(i_l).wXww;
        avU(i_l).wXS2  = sqrt(avU(i_l).wXuu.^2 + avU(i_l).wXvv.^2 + avU(i_l).wXww.^2);
    end
    
    % Consistency check for wB = wBt + wBq
    if all(ismember({'wBt','wBq'},fieldnames(avS)))
        avS(i_l).wQ  = avS(i_l).wBt + avS(i_l).wBq;
        avU(i_l).wQ  = sqrt(avU(i_l).wBt.^2 + avU(i_l).wBq.^2);
    end
    
end



%% Compensate / integrate

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    r = (1:length(avS(i_l).uuu))*avS(i_l).dr;
    
    avS(i_l).W  = 6*cumtrapz(r,avS(i_l).wB .*r.^2) ./ r.^3;
    avS(i_l).Wt = 6*cumtrapz(r,avS(i_l).wBt.*r.^2) ./ r.^3;
    avS(i_l).Wq = 6*cumtrapz(r,avS(i_l).wBq.*r.^2) ./ r.^3;
    avU(i_l).W  = 6*cumtrapz(r,avU(i_l).wB .*r.^2) ./ r.^3;
    avU(i_l).Wt = 6*cumtrapz(r,avU(i_l).wBt.*r.^2) ./ r.^3;
    avU(i_l).Wq = 6*cumtrapz(r,avU(i_l).wBq.*r.^2) ./ r.^3;
    
    avS(i_l).Ti = 3*cumtrapz(r,avS(i_l).wXS2.*r.^2) ./ r.^3;
    avU(i_l).Ti = 3*cumtrapz(r,avU(i_l).wXS2.*r.^2) ./ r.^3;
    
    avS(i_l).uuu3r = 3*avS(i_l).uuu./r;
    avS(i_l).vvu3r = 3*avS(i_l).vvu./r;
    avS(i_l).wwu3r = 3*avS(i_l).wwu./r;
    avU(i_l).uuu3r = 3*avU(i_l).uuu./r;
    avU(i_l).vvu3r = 3*avU(i_l).vvu./r;
    avU(i_l).wwu3r = 3*avU(i_l).wwu./r;
    
    avS(i_l).S3r = avS(i_l).S3./r;
    avU(i_l).S3r = avU(i_l).S3./r;
    
    avS(i_l).mS3 = -avS(i_l).S3;
    avU(i_l).mS3 = avU(i_l).S3;
    
    avS(i_l).WmS3r = avS(i_l).W - avS(i_l).S3r;
    avU(i_l).WmS3r = avU(i_l).W + avU(i_l).S3r;
    
    avS(i_l).WrmS3 = avS(i_l).W.*r - avS(i_l).S3;
    avU(i_l).WrmS3 = avU(i_l).W.*r + avU(i_l).S3;
        
end



%% Estimate dissipation rate

edr_vars = {'uu','vv','ww','mS3','WrmS3'};
slps = [2/3 2/3 2/3 1 1 1];
cons = [2.0 2.6 2.6 4.0 4.0 4.0];

Nvar = numel(edr_vars);

for i_v = 1:Nvar
    var = edr_vars{i_v};
    [avS,avU] = calc_edr(avS,avU,edr_vars{i_v},edr_fit_range,cons(i_v),...
        'Slope',slps(i_v),'Method','direct');
end


% Compensate for plotting

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    r = (1:length(avS(i_l).uuu))*avS(i_l).dr;
    
    for i_v = 1:Nvar
        var = edr_vars{i_v};
        avS(i_l).([var,'_c']) = avS(i_l).(var) ./ r.^slps(i_v) / cons(i_v);
        avU(i_l).([var,'_c']) = avU(i_l).(var) ./ r.^slps(i_v) / cons(i_v);
    end
    
    avS(i_l).edr_S2_4 = avS(i_l).edr_S2*4;
    avU(i_l).edr_S2_4 = avU(i_l).edr_S2*4;
end



%% Save/load

save('S_baltic.mat','r_max','edr_fit_range','levels',...
    'MOM','S','U','L','N','V','avS','avU','avN',...
    'plotpath')


% load('S_baltic.mat')



%% TABLES

% Segment overview

print_table(MOM,{'length_km','alt'},true,0,{'mean','std'})


% Dissipation rates

print_vars = {'uu','vv','ww','S2','mS3','WrmS3'};%,'WrmS3mTr'};

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

xlim = [dr r_max];
ylim = [1e-8 1e-1];

Npoints = 40;

Nlvl = size(avS,1);


% S2

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'uu_c','vv_c','ww_c'},[7 5 1],Npoints,'XLim',xlim);
    if i_l>0
        legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    end
    ylabel('$S_2r^{-2/3}C^{-1} $','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s2_',levels{i_l}],'-dpng','-r300')
end
    

% W

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','Wt','Wq'},[7 5 1],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$W$','$W_\theta$','$W_q$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'W_',levels{i_l}],'-dpng','-r300')
end


% S3

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'S3r'},[1 6],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$S_3 r^{-1}$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s3r_',levels{i_l}],'-dpng','-r300')
end


% Balance: W, S3/r, epsilon

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','-S3r','edr_S2_4'},[7 1 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$W$','$-S_3 r^{-1}$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'balance_',levels{i_l}],'-dpng','-r300')
end

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','-S3r','edr_S2_4'},[7 1 3],Npoints,'XLim',xlim,'YScale','linear');
    if i_l>0
        legend({'$W$','$-S_3 r^{-1}$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'balance_lin_',levels{i_l}],'-dpng','-r300')
end


% Total sum: W - S3/r

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'WmS3r','edr_S2_4'},[4 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$W-S_3 r^{-1}$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt)) 
    print(fig,[plotpath,filesep,'total_',levels{i_l}],'-dpng','-r300')
end

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'WmS3r','edr_S2_4'},[4 3],Npoints,'XLim',xlim,'YScale','linear');
    if i_l>0
        legend({'$W-S_3 r^{-1}$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt)) 
    print(fig,[plotpath,filesep,'total_lin_',levels{i_l}],'-dpng','-r300')
end


% S3 components

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'uuu3r','vvu3r','wwu3r'},[7 5 1 2 4 6],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'uuu','vvu','wwu'},'Interpreter','latex','Location','best')
    end
    ylabel('$3S_3 r^{-1}\,[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s3_comp_',levels{i_l}],'-dpng','-r300')
end


% Fitting dissipation

for i_l = 1:Nlvl
    fig = plot_sfc_edr(avS(i_l),[],{'uu','vv','ww'},[7 5 1],Npoints,'XLim',xlim);
    if i_l>0
        legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    end
    ylabel('$S_2r^{-2/3}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr2fit_',levels{i_l}],'-dpng','-r300')
end

for i_l = 1:Nlvl
    fig = plot_sfc_edr(avS(i_l),[],{'mS3','WrmS3'},[1 7 4],Npoints,'XLim',xlim);
    if i_l>0
        legend({'$-S_3 r^{-1}$','$W-S_3 r^{-1}$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^3\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr3fit_',levels{i_l}],'-dpng','-r300')
end