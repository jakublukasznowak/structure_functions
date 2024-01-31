
% MYPROJECTPATH is the path where you downloaded the codes
%
% MYDATAPATH is the path where you downloaded two datasets:
%
% MYDATAPATH/TURBLENCE
%
% Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz. Aeris.
% https://doi.org/10.25326/128
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/
% 
% In this code 'longlegs' L3 v1.9 is used.
% 
% MYDATAPATH/CloudComposite
%
% Coutris, P. (2021). SAFIRE ATR42: PMA/Cloud composite dataset. Aeris.
% https://doi.org/10.25326/237
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/PMA/PROCESSED/CloudComposite/
% 
% In this code v1 is used.


% Prepare paths

addpath(genpath(myprojectpath))

datapath = mydatapath;

plotpath = [myprojectpath,filesep,'figures'];
if ~isfolder(plotpath)
    mkdir(plotpath)
end



%% Load datasets


% List of levels
levels  = {'cloud-base','cloud-base-noclouds','top-subcloud','mid-subcloud','near-surface'};

% List of flights
flight_ids = num2cell(num2str((9:19)','RF%02d'),2); % RF09 - RF19

% List of variables from turbulent moments dataset
mom_vars = {'alt','time_start','time_end',...
    'MEAN_P','MEAN_THETA','MEAN_MR',...
    'MEAN_WSPD','MEAN_WDIR','MEAN_TAS','MEAN_THDG'};

% List of variables from turbulent fluctuations dataset
turb_vars = {'time','W_DET','UX_DET','VY_DET','T_DET','MR_DET'};

% List of variables from cloud composite dataset
pma_vars = {'time','LWC','NT','CLOUD_mask'};


% Read data files

% Flight segmentation
disp('Load flight segmentation ...')
SEG = load_atr_seg(datapath,'v1.9','longlegs');                          
SEG = SEG(ismember(SEG.flight,flight_ids) & ismember(SEG.level,levels),:);

% Mean values and moments
disp('Load mean values and moments ...')
[MOM,mom_info] = load_atr_mom(datapath,'L3','v1.9','longlegs',mom_vars); 
MOM = join(SEG,MOM,'Keys',{'start','end'});

% Turbulent fluctuations
disp('Load turbulence data ...')
[TURB,turb_info] = load_atr_turb(MOM,datapath,'L3','v1.9',turb_vars);

% Microphysics
disp('Load microphysics data ...')
PMA = load_atr_pma(MOM,datapath,pma_vars);


% Calculate additional parameters

% Thermodynamics (including buoyancy)
TURB = calc_thermo(TURB,MOM);

% Cloud masks
% (1) from PMA dataset based on LWC and interpolated to 25 Hz
% (2) based on RH>98%
% [both extended in front and behind each cloud by its 1 width]
% (3) union of the two extended masks
[TURB,MOM] = cloud_mask(TURB,MOM,PMA,98,1);

MOM.valid_fraction = 1 - MOM.OR_cloud_fraction;



%% Calculate SFC with uncertainties for individual segments

% List of structure functions
sfc_vars = {'uuu',{'VY_DET','VY_DET','VY_DET'};
            'vvu',{'UX_DET','UX_DET','VY_DET'};
            'wwu',{'W_DET','W_DET','VY_DET'};
            'wB', {'W_DET','B'};
            'wBt',{'W_DET','Bt'};
            'wBq',{'W_DET','Bq'};
            'uu', {'VY_DET','VY_DET'};
            'vv', {'UX_DET','UX_DET'};
            'ww', {'W_DET','W_DET'};
            'uuW',{'VY_DET','VY_DET','X_W_DET'};
            'vvW',{'UX_DET','UX_DET','X_W_DET'};
            'wwW',{'W_DET','W_DET','X_W_DET'}
            };

% Assume fixed distance between points
MOM.dr(:) = 4;

% Consider max lag of 400 m
r_maxlag = 400/4;

% Duplicate cloud-base segments
Nseg = size(MOM,1);
MOM2 = [MOM; MOM(MOM.level=="cloud-base",:)];
MOM2.level(Nseg+1:end) = "cloud-base-noclouds";
TURB2 = [TURB; TURB([TURB(:).level]=="cloud-base")];
MOM2.valid_fraction(MOM.level~="cloud-base-noclouds") = 1;


S = table2struct(MOM2(:,{'flight','name','level','alt','dr','length','valid_fraction'}));
V = S; L = S; N = S; U = S;

Nseg = size(S,1);
Nvar = size(sfc_vars,1);

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
        
        Lt = length(TURB2(i_s).time);
    
        A = cell2mat(cellfun(@(v) TURB2(i_s).(v),varlist,'UniformOutput',false));
        
        % Mask cloudy points
        if S(i_s).level == "cloud-base-noclouds"
            A(TURB2(i_s).OR_mask,:) = nan;
        end
        
        if ~any(all(isnan(A),1))
            
            S(i_s).(var) = nan(1,r_maxlag);
            V(i_s).(var) = nan(1,r_maxlag);
            L(i_s).(var) = nan(1,r_maxlag);
            S(i_s).N     = nan(1,r_maxlag);

            for r_lag = 1:r_maxlag
                
                difI = prod( A(r_lag+1:Lt,difFlag)-A(1:Lt-r_lag,difFlag) ,2);
                
                if sum(midFlag)>0
                    midI = prod( 0.5*(A(r_lag+1:Lt,midFlag)+A(1:Lt-r_lag,midFlag)) ,2);
                else
                    midI = 1;
                end
                
                I = difI.*midI; % increments
                
                S(i_s).(var)(r_lag) = mean(I,'omitnan');        % structure function
                V(i_s).(var)(r_lag) = mean(I.^2,'omitnan');     % variance of increments
                L(i_s).(var)(r_lag) = int_ls_short(I)*S(i_s).dr;% integral length scale for increments
                S(i_s).N    (r_lag) = sum(~isnan(I));           % number of valid samples
            end
            
            % number of independent samples
%             N(i_s).(var) = S(i_s).N*S(i_s).dr./L(i_s).(var); 
            N(i_s).(var) = S(i_s).valid_fraction*S(i_s).length./L(i_s).(var); 
            
            % uncertainty
            U(i_s).(var) = sqrt( (V(i_s).(var) - S(i_s).(var).^2) ./ N(i_s).(var) ); 
              
        end

    end
    
end



%% Average SFC at 4 characteristic levels

% List of segments to exclude from averaging
% chosen based on the manual inspection of the timeseries of altitude,
% temperature, wind speed and direction, LWC, TAS, heading from 
% the CORE dataset (not shown here):
% 
% CNRM/TRAMM, SAFIRE, Laboratoire d'Aérologie. (2021). SAFIRE ATR42: Core 
% Data 1Hz - V2. Aeris. https://doi.org/10.25326/298

exclude_seg = ["RF19","S1";   % (surface) altitude changes + temperature drop
               "RF09","L1B";  % (mid-subcloud) temperature drop + wind direction variations + TAS changes
               "RF11","L2B";  % (mid-subcloud) rain + TAS changes NEW
               "RF17","L2A";  % (mid-subcloud) intensive rain + temperature drop + TAS changes
               "RF11","L1A";  % (top-subcloud) intensive rain
               "RF17","L1B";  % (top-subcloud) intensive rain + temperature variations + TAS changes
               "RF17","R2B";  % (cloud-base) heavy rain NEW
               "RF18","R1B"]; % (cloud-base) heavy rain NEW
 

Nseg = size(S,1);
Nlvl = numel(levels);
Nvar = size(sfc_vars,1);

% Valid segments
valid_seg = true(Nseg,1);
for i_e = 1:size(exclude_seg,1)
    valid_seg( [S(:).flight]==exclude_seg(i_e,1) & [S(:).name]==exclude_seg(i_e,2) ) = false;
end

% Segments with correct humidity measurement
valid_hum = (MOM2.RH_nan_fraction==0);


avS = struct([]); avN = struct([]); avU = struct([]);

for i_l = 1:Nlvl
    ind_l = find([S(:).level]'==levels{i_l} & valid_seg);
    
    avS(i_l,1).level  = string(levels{i_l});
    avS(i_l).number = numel(ind_l);
    avS(i_l).alt    = mean([S(ind_l).alt]);
    avS(i_l).dr     = mean([S(ind_l).dr]);
    avS(i_l).length = mean([S(ind_l).length]);
    
    for i_v = 1:Nvar
        var = sfc_vars{i_v,1};
        
        if ismember(var,{'wB','wBt','wBq'})
            ind_lv = find([S(:).level]'==levels{i_l} & valid_seg & valid_hum);
        else
            ind_lv = ind_l;
        end
        
        avS(i_l,1).(var) = mean( cat(1,S(ind_lv).(var)), 1);
        avV              = mean( cat(1,V(ind_lv).(var)), 1);
%         avN(i_l,1).(var) = sum( cat(1,N(ind_lv).(var)) ./ cat(1,L(ind_lv).(var)) );
        avN(i_l,1).(var) = sum( repmat([S(ind_lv).valid_fraction]'.*[S(ind_lv).length]',1,r_maxlag) ...
            ./ cat(1,L(ind_lv).(var)) );
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
        avU(i_l).s2W  = sqrt(avS(i_l).uuW.^2 + avS(i_l).vvW.^2 + avS(i_l).wwW.^2);
    end
end



%% Compensate / integrate

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    r = (1:r_maxlag)*avS(i_l).dr;
    
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
    
    avS(i_l).Ti = 3*cumtrapz(r,avS(i_l).s2W.*r.^2) ./ r.^3;
    avU(i_l).Ti = 3*cumtrapz(r,avU(i_l).s2W.*r.^2) ./ r.^3;
    
%     avS(i_l).uuu3r = 3*avS(i_l).uuu./r;
%     avS(i_l).vvu3r = 3*avS(i_l).vvu./r;
%     avS(i_l).wwu3r = 3*avS(i_l).wwu./r;
%     avU(i_l).uuu3r = 3*avU(i_l).uuu./r;
%     avU(i_l).vvu3r = 3*avU(i_l).vvu./r;
%     avU(i_l).wwu3r = 3*avU(i_l).wwu./r;
    
    avS(i_l).Wrs3l = avS(i_l).W.*r - avS(i_l).s3l;
    avU(i_l).Wrs3l = avU(i_l).W.*r + avU(i_l).s3l;
end



%% Calculate dissipation rate

fit_range = [20 60];

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



%% Estimate transport

% as residual

Nlvl = size(avS,1);
for i_l = 1:Nlvl
    avS(i_l).Tres = -avS(i_l).s3lr + avS(i_l).W - 4*avS(i_l).edr_s2;
    avU(i_l).Tres = sqrt( avU(i_l).s3lr.^2 + avU(i_l).W.^2 + 4*avU(i_l).edr_s2.^2 );
end

% from differences between levels

level_pairs = {'cloud-base','top-subcloud';
               'cloud-base-noclouds','top-subcloud';
               'top-subcloud','mid-subcloud';
               'mid-subcloud','near-surface'};

for i_l = 1:size(level_pairs,1)
    i1 = find([avS(:).level]==level_pairs{i_l,1});
    i2 = find([avS(:).level]==level_pairs{i_l,2});
    
    avS(i1).Tdif = (avS(i1).Ti-avS(i2).Ti) / (avS(i1).alt-avS(i2).alt);
    avU(i1).Tdif = sqrt(avU(i1).Ti.^2+avU(i2).Ti.^2) / (avS(i1).alt-avS(i2).alt);
end



%% TABLES

% Segment overview

print_table(MOM2(valid_seg,:),{'length_km','alt'},1,0)


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

% save('S.mat','S','U','N','L','V','avS','avU','avN','avN','avV')
% load('S.mat')


%% Overview of the segments

plot_seg_overview(MOM,setdiff(levels,{'cloud-base-noclouds'},'stable'),false);
print(gcf,[plotpath,filesep,'seg_overview'],'-dpng','-r300')


%% Structure functions

plots = {
%     {'uuu3r','vvu3r','wwu3r'}, {'uuu','vvu','wwu'}, '$[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
    {'W','Wt','Wq'}, {'$W$','$W_\theta$','$W_q$'},'$W\,[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
    {'s3lr'}, {'$S_3^L r^{-1}$'}, '$S_3^Lr^{-1}\,[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
    {'Ws3lr','edr_s2_4'}, {'$W-S_3^L r^{-1}$','$4\epsilon_2$'},'$[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
    {'Tres','Tdif'}, {'$T_{res}$','$T_{dif}$'}, '$T\,[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
    {'uu_c','vv_c','ww_c'}, {'uu','vv','ww'}, '$S_2r^{-2/3}C^{-1}$', [8e-4 1e-2];
    };


Nlvl = size(avS,1);
      
for i_p = 1:size(plots,1)
    for i_l = 1:Nlvl
        fig = plot_sfc(avS(i_l),avU(i_l),plots{i_p,1},[7 1 5],...
            'XLim',[4 400],'YLim',plots{i_p,4});
        if numel(plots{i_p,2})>1
            legend(plots{i_p,2},'Interpreter','latex','Location','best')
        end
        ylabel(plots{i_p,3},'Interpreter','latex')
        title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
        print(fig,[plotpath,filesep,plots{i_p,1}{1},'_',levels{i_l}],'-dpng','-r300')
    end
end



%% Fitting dissipation

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    
    fig = plot_sfc_edr(avS(i_l),[],{'uu','vv','ww'},[5 6 7],'XLim',[4 400],'YLim',[1e-3 2e-2]);
    legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    ylabel('$S_2r^{-2/3}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr2_',levels{i_l}],'-dpng','-r300')
    
    fig = plot_sfc_edr(avS(i_l),avU(i_l),{'Wrs3l'},7,'XLim',[4 400],'YLim',[1e-6 2e-3]);
    ylabel('$W-S_3^L r^{-1}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr3_',levels{i_l}],'-dpng','-r300')
    
end



%% OLD

% Stare obliczenia wX - niepoprawne
%
% Calculate velocity in between data points W_MID 
% Nseg = size(TURB,1);
% for i_s = 1:Nseg
%     Lt = length(TURB(i_s).time);
%     TURB(i_s).W_MID = interp1(TURB(i_s).W_DET,0.5:0.5:Lt)';
%     TURB(i_s).OR_mask_mid = (interp1(double(TURB(i_s).OR_mask),0.5:0.5:Lt)'>0);
% end
% (...)        
%         La = cellfun(@length,A);
%         difA = cell2mat(A(La==Lt));
%         midA = cell2mat(A(La>Lt));
%         if S(i_s).level == "cloud-base"
%             difA(TURB(i_s).OR_mask,:) = nan; 
%             if ~isempty(midA)
%                 midA(TURB(i_s).OR_mask_mid,:) = nan; 
%             end
%         end
% (...)
%          midI = prod( midA(r_lag+1:2:2*Lt-r_lag,:) ,2);
