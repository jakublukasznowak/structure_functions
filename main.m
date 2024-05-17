
%% Introduction

% The code requires access to the datasets which can be downloaded from
% the public data repositories.
%
% MYPROJECTPATH is the location where you downloaded the codes
% MYDATAPATH is the path where you downloaded the datasets
%
%
% MYDATAPATH/ATR-EUREC4A/TURBLENCE
%
% Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz. Aeris.
% doi.org/10.25326/128
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/
% 
% In this code 'longlegs' L3 v1.9 is used.
% 
%
% MYDATAPATH/CloudComposite
%
% Coutris, P. (2021). SAFIRE ATR42: PMA/Cloud composite dataset. Aeris.
% https://doi.org/10.25326/237
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/PMA/PROCESSED/CloudComposite/
% 
% In this code v1 is used.
%
%
%
% The code was developed in MATLAB R2019b. The functionality in other
% versions of this environment was not tested.
%
% One external packages is used:
%       (1) YAML 1.1 parser and emitter for MATLAB by Martin Koch
%           from https://www.mathworks.com/matlabcentral/fileexchange/106765-yaml



%% Settings

% Compute uncertainties?
% (takes quite a long time)

if_uncertainties = false;


% Cloud mask

mask_rh_threshold = 98;   % Relative humidity threshold for cloud mask
mask_ext_factor   = 1;    % Mask extension factor to account for cloud shells
mask_nan_replace  = true; % Assume cloud if masks are nan


% Structure functions

r_max = 1000; % [m] Max r separation


% Transport term

transport_method = 'interp'; % interp or diff


% Dissipation rate

edr_fit_range = [20 60]; % [m]



%% Load list

% List of levels
levels  = {'cloud-base','cloud-base-noclouds','top-subcloud','mid-subcloud','near-surface'};

% List of flights
flight_ids = num2cell(num2str((9:19)','RF%02d'),2); % RF09 - RF19

% List of variables from turbulent moments dataset
mom_vars = {'alt','time_start','time_end',...
    'MEAN_P','MEAN_THETA','MEAN_MR','MEAN_WDIR','MEAN_TAS','MEAN_THDG'};

% List of variables from turbulent fluctuations dataset
turb_vars = {'time','UX_DET','VY_DET','W_DET','T_DET','MR_DET'};

% List of variables from cloud microphysics dataset
mcph_vars = {'time','LWC','CLOUD_mask'};


% List of segments to exclude
% (chosen based on the manual inspection of the timeseries of altitude,
% temperature, wind speed and direction, LWC, TAS and true heading
% from the CORE dataset:
% CNRM/TRAMM, SAFIRE, Laboratoire d'Aérologie. (2021). SAFIRE ATR42: Core 
% Data 1Hz - V2. Aeris. https://doi.org/10.25326/298

exclude_seg = ["RF19","S1";   % (surface) altitude changes + temperature drop
               "RF09","L1B";  % (mid-subcloud) temperature drop + wind direction variations + TAS changes
               "RF11","L2B";  % (mid-subcloud) rain + TAS changes
               "RF17","L2A";  % (mid-subcloud) heavy rain + temperature drop + TAS changes
               "RF11","L1A";  % (top-subcloud) heavy rain
               "RF17","L1B";  % (top-subcloud) heavy rain + temperature variations + TAS changes
               "RF17","R2B";  % (cloud-base) heavy rain
               "RF18","R1B"]; % (cloud-base) heavy rain

           

%% Prepare paths

addpath(genpath(myprojectpath))

datapath = mydatapath;

plotpath = [myprojectpath,filesep,'figures'];
if ~isfolder(plotpath), mkdir(plotpath), end



%% Load datasets

% Flight segmentation

disp('Load flight segmentation ...')
SEG = load_atr_seg(datapath,'v1.9','longlegs');                          
SEG = SEG(ismember(SEG.flight,flight_ids) & ismember(SEG.level,levels),:);


% Exclude problematic segments

Nexc = size(exclude_seg,1);
for i_e = 1:Nexc
    ind_s = find( SEG.flight==exclude_seg(i_e,1) & SEG.name==exclude_seg(i_e,2) );
    SEG(ind_s,:) = [];
end


% Mean values and moments

disp('Load mean values and moments ...')
[MOM,mom_info] = load_atr_mom(datapath,'L3','v1.9','longlegs',mom_vars); 
MOM = join(SEG,MOM,'Keys',{'start','end'});


% Turbulent fluctuations

disp('Load turbulence data ...')
[TURB,turb_info] = load_atr_turb(MOM,datapath,'L3','v1.9',turb_vars);

MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';


% Microphysics

disp('Load microphysics data ...')
MCPH = load_atr_pma(MOM,datapath,mcph_vars);


Nseg = size(MOM,1);
fprintf('Number of segments loaded: %d\n',Nseg)



%% Compute auxilliary parameters

disp('Compute thermodynamics and cloud mask ...')


% Thermodynamics (including buoyancy)

TURB = calc_thermo(TURB,MOM);


% Cloud masks

% (1) from microphysics based on LWC and interpolated to 25 Hz
% (2) based on RH>threshold
% [both extended in front and behind each cloud]
% (3) union of the two extended masks

[TURB,MOM] = cloud_mask(TURB,MOM,MCPH,mask_rh_threshold,mask_ext_factor,mask_nan_replace);

for i_s = 1:Nseg
    TURB(i_s).OR_mask(isnan(TURB(i_s).OR_mask)) = mask_nan_replace;
end
MOM.valid_fraction = 1 - MOM.OR_cloud_fraction;


% Duplicate cloud-base segments

MOM = [MOM; MOM(MOM.level=="cloud-base",:)];
MOM.level(Nseg+1:end) = "cloud-base-noclouds";
MOM.valid_fraction(MOM.level~="cloud-base-noclouds") = 1;

TURB = [TURB; TURB([TURB(:).level]=="cloud-base")];


% Assume fixed distance between points

dr = 4;
MOM.dr(:) = dr;
r_maxlag = r_max/dr;



%% Compute structure functions for individual segments

disp('Compute structure functions ...')

% List of structure functions
sfc_vars = {'uuu',{'UX_DET','UX_DET','UX_DET'};
            'vvu',{'VY_DET','VY_DET','UX_DET'};
            'wwu',{'W_DET','W_DET','UX_DET'};
            'wB', {'W_DET','B'};
            'wBt',{'W_DET','Bt'};
            'wBq',{'W_DET','Bq'};
            'uu', {'UX_DET','UX_DET'};
            'vv', {'VY_DET','VY_DET'};
            'ww', {'W_DET','W_DET'};
            'wXuu',{'UX_DET','UX_DET','X_W_DET'};
            'wXvv',{'VY_DET','VY_DET','X_W_DET'};
            'wXww',{'W_DET','W_DET','X_W_DET'}};


% Prepare data structures

S = table2struct(MOM(:,{'flight','name','level','alt','dr','length','valid_fraction'})); % structure functions
V = S; % Variances of the form <dx^2>
L = S; % Integral length scales of increments
N = S; % Number of independent samples
U = S; % Uncertainties of structure functions


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
        
        % Data matrix
        A = cell2mat(cellfun(@(v) TURB(i_s).(v),varlist,'UniformOutput',false));
        
        % Mask cloudy points
        if S(i_s).level == "cloud-base-noclouds"
            A(logical(TURB(i_s).OR_mask),:) = nan;
        end
        
        Lt = length(TURB(i_s).time);
        
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
%             N(i_s).(var) = S(i_s).N*S(i_s).dr./L(i_s).(var); 
            N(i_s).(var) = S(i_s).valid_fraction*S(i_s).length./L(i_s).(var); 
            
            % uncertainty
            U(i_s).(var) = sqrt( (V(i_s).(var) - S(i_s).(var).^2) ./ N(i_s).(var) ); 
              
        end
        
        us = etd(us,Nseg*(i_v-1)+i_s);

    end
    
end



%% Average at characteristic levels

disp('Average at characteristic levels ...')


% Segments with correct humidity measurement
valid_hum = (MOM.RH_nan_fraction==0);


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
        
        if ismember(var,{'wB','wBt','wBq'})
            ind_lv = find([S(:).level]'==levels{i_l} & valid_hum);
        else
            ind_lv = ind_l;
        end
        
        avS(i_l,1).(var) = mean( vertcat(S(ind_lv).(var)), 1);
        avV              = mean( vertcat(V(ind_lv).(var)), 1);
%         avN(i_l,1).(var) = sum( vertcat(N(ind_lv).(var)) ./ vertcat(L(ind_lv).(var)) );
        avN(i_l,1).(var) = sum( repmat([S(ind_lv).valid_fraction]'.*[S(ind_lv).length]',1,r_maxlag) ...
            ./ vertcat(L(ind_lv).(var)) );
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



%% Integrate / compensate

disp('Integrate or compensate ...')

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



%% Estimate transport

disp('Estimate transport ...')


if strcmp(transport_method,'diff') % from differences between levels
    
    level_pairs = {'cloud-base','top-subcloud';
                   'cloud-base-noclouds','top-subcloud';
                   'top-subcloud','mid-subcloud';
                   'mid-subcloud','near-surface'};

    for i_lp = 1:size(level_pairs,1)
        i_l1 = find([avS(:).level]==level_pairs{i_lp,1});
        i_l2 = find([avS(:).level]==level_pairs{i_lp,2});

        avS(i_l1).T = (avS(i_l1).Ti-avS(i_l2).Ti) / (avS(i_l1).alt-avS(i_l2).alt);
        avU(i_l1).T = sqrt(avU(i_l1).Ti.^2+avU(i_l2).Ti.^2) / (avS(i_l1).alt-avS(i_l2).alt);
    end

    
elseif strcmp(transport_method,'interp') % from interpolated smooth function

    alt_increment = 25; % m
    interp_method = 'pchip';

    level_sets = {'cloud-base','top-subcloud','mid-subcloud','near-surface';
         'cloud-base-noclouds','top-subcloud','mid-subcloud','near-surface'};
    
    for i_ls = 1:size(level_sets,1) 
        [~,ind_l] = ismember(level_sets(i_ls,:),[avS(:).level]);
        
        alt_vec   = [avS(ind_l).alt];
        alt_query = reshape( [alt_vec-alt_increment; alt_vec+alt_increment] ,1,[]);
        
        Ti_interp = interp1(alt_vec,vertcat(avS(ind_l).Ti),alt_query,interp_method);
        if if_uncertainties
            u_Ti_interp = interp1(alt_vec,cat(1,avU(ind_l).Ti),alt_query,interp_method);
        else
            u_Ti_interp = nan(size(Ti_interp));
        end

        for ii_l = 1:numel(ind_l)
            i_l = ind_l(ii_l);
            avS(i_l).T = ( Ti_interp(2*ii_l,:) - Ti_interp(2*ii_l-1,:) )/alt_increment/2;
            avU(i_l).T = sqrt( u_Ti_interp(2*ii_l,:).^2 + u_Ti_interp(2*ii_l-1,:).^2 )/alt_increment/2;
        end
    end
    
end


% Dependent parameters

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    r = (1:length(avS(i_l).uuu))*avS(i_l).dr;

    avS(i_l).WmS3rmT = avS(i_l).W - avS(i_l).S3r - avS(i_l).T;
    avU(i_l).WmS3rmT = avU(i_l).W + avU(i_l).S3r + avU(i_l).T;
    
    avS(i_l).WrmS3mTr = avS(i_l).W.*r - avS(i_l).S3 - avS(i_l).T.*r;
    avU(i_l).WrmS3mTr = avU(i_l).W.*r + avU(i_l).S3 + avU(i_l).T.*r;
end



%% Estimate dissipation rate

disp('Estimate dissipation rate ...')

edr_vars = {'uu','vv','ww','mS3','WrmS3','WrmS3mTr'};
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

% save('S_eureca.mat','r_max','transport_method','edr_fit_range','levels','exclude_seg','mask*',...
%     'MOM','S','U','L','N','V','avS','avU','avN',...
%     'plotpath')


% load('S_eureca.mat')

          

%% TABLES

% Segment overview

print_table(MOM,{'length_km','alt'},true,0,{'mean','std'})


% Dissipation rates

print_vars = {'uu','vv','ww','S2'};%'mS3','WrmS3','WrmS3mTr'};

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

xlim = [4 1000];
ylim = [1e-6 5e-3];

Npoints = 40;

Nlvl = size(avS,1);


% S2

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'uu_c','vv_c','ww_c'},[7 5 1],Npoints,'XLim',xlim,'YLim',[8e-4 1e-2]);
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


% Balance: W, S3/r, T, epsilon

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','-S3r','-T','edr_S2_4'},[7 1 5 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$W$','$-S_3 r^{-1}$','$-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'balance_',levels{i_l}],'-dpng','-r300')
end

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','-S3r','-T','edr_S2_4'},[7 1 5 3],Npoints,'XLim',xlim,'YScale','linear');
    if i_l>0
        legend({'$W$','$-S_3 r^{-1}$','$-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'balance_lin_',levels{i_l}],'-dpng','-r300')
end


% Total sum: W - S3/r - T

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'WmS3rmT','edr_S2_4'},[4 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l>0
        legend({'$W-S_3 r^{-1}-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt)) 
    print(fig,[plotpath,filesep,'total_',levels{i_l}],'-dpng','-r300')
end

for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'WmS3rmT','edr_S2_4'},[4 3],Npoints,'XLim',xlim,'YScale','linear');
    if i_l>0
        legend({'$W-S_3 r^{-1}-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
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
    fig = plot_sfc_edr(avS(i_l),[],{'uu','vv','ww'},[7 5 1],Npoints,'XLim',xlim,'YLim',[1e-3 2e-2]);
    if i_l>0
        legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    end
    ylabel('$S_2r^{-2/3}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr2fit_',levels{i_l}],'-dpng','-r300')
end

for i_l = 1:Nlvl
    fig = plot_sfc_edr(avS(i_l),[],{'mS3','WrmS3','WrmS3mTr'},[1 7 4],Npoints,'XLim',xlim,'YLim',[1e-5 1e-2]);
    if i_l>0
        legend({'$-S_3 r^{-1}$','$W-S_3 r^{-1}$','$W-S_3 r^{-1}-T_u$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^3\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr3fit_',levels{i_l}],'-dpng','-r300')
end


% Check power law in T

t_levels = {'cloud-base','cloud-base-noclouds','top-subcloud'};
t_fit_range = [8 80; 8 300; 8 100];


for ii_l = 1:numel(t_levels)
    i_l = find(strcmp(levels,t_levels{ii_l}));
    
    y = -avS(i_l).T;
    r = (1:length(y))*S(i_l).dr;

    ind_r = unique(interp1(r,1:length(y),exp(linspace(log(t_fit_range(i_l,1)),log(t_fit_range(i_l,2)),Npoints)),'nearest',length(y)));
    rfit = r(ind_r);
    yfit = y(ind_r);
    
    [p,SM] = polyfit(log(rfit),log(yfit),1);

    slp = p(1);
    logOstar = p(2);
    
    covarM = (inv(SM.R)*inv(SM.R)')*SM.normr^2/SM.df; % covariance matix
    e.slp  = sqrt(covarM(1,1));
    e.logOstar = sqrt(covarM(2,2));
    
    Ostar = exp(logOstar);
    e.Ostar = Ostar*e.logOstar;
    
    corrM = corrcoef(log(rfit),log(yfit)); % correlation matix
    e.R2 = corrM(1,2);

    [fig,ax] = plot_sfc(avS(i_l),avU(i_l),{'-T'},5,Npoints,'XLim',xlim);%,'YLim',ylim);
    plot(ax,rfit,Ostar*rfit.^slp,'-','Color','blue','LineWidth',2)
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    legend({'$-T_u$',['$s=$ ',num2str(slp,'%.2f')]},'Interpreter','latex','Location','best')
    print(fig,[plotpath,filesep,'Tfit_',levels{i_l}],'-dpng','-r300')
end