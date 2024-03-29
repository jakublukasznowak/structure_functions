
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

iferrors = true;
orientation = 'reversed';


% Prepare paths

addpath(genpath(myprojectpath))

datapath = mydatapath;

plotpath = [myprojectpath,filesep,'figures_',orientation];
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

if strcmp(orientation,'reversed')
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
                'wwW',{'W_DET','W_DET','X_W_DET'};
%                 'uuA',{'VY_DET','VY_DET','A_VY_DET'};
%                 'vvA',{'UX_DET','UX_DET','A_VY_DET'};
%                 'wwA',{'W_DET','W_DET','A_VY_DET'};
                };
elseif strcmp(orientation,'original')
    sfc_vars = {'uuu',{'UX_DET','UX_DET','UX_DET'};
                'vvu',{'VY_DET','VY_DET','UX_DET'};
                'wwu',{'W_DET','W_DET','UX_DET'};
                'wB', {'W_DET','B'};
                'wBt',{'W_DET','Bt'};
                'wBq',{'W_DET','Bq'};
                'uu', {'UX_DET','UX_DET'};
                'vv', {'VY_DET','VY_DET'};
                'ww', {'W_DET','W_DET'};
                'uuW',{'UX_DET','UX_DET','X_W_DET'};
                'vvW',{'VY_DET','VY_DET','X_W_DET'};
                'wwW',{'W_DET','W_DET','X_W_DET'};
%                 'uuA',{'UX_DET','UX_DET','A_UX_DET'};
%                 'vvA',{'VY_DET','VY_DET','A_UX_DET'};
%                 'wwA',{'W_DET','W_DET','A_UX_DET'};
                };
end
            
% Max r separation
r_max = 1000;
dr = 4;

% Assume fixed distance between points
MOM.dr(:) = dr;
r_maxlag = r_max/4;

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
us = etd(clock,0,Nseg*Nvar,30);

for i_v = 1:Nvar
    
    var = sfc_vars{i_v,1};
    varlist = sfc_vars{i_v,2};
    
    midFlag = startsWith(varlist,'X_');
    absFlag = startsWith(varlist,'A_');
    difFlag = ~startsWith(varlist,'X_') & ~startsWith(varlist,'A_');
    varlist(midFlag) = extractAfter(varlist(midFlag),'X_');
    varlist(absFlag) = extractAfter(varlist(absFlag),'A_');
    
    fprintf('S_%s = < ',var)
    if sum(midFlag)
        fprintf('%s_X ',varlist{midFlag})
    end
    if sum(absFlag)
        fprintf('|d %s| ',varlist{absFlag})
    end
    if sum(difFlag)
        fprintf('d %s ',varlist{difFlag})
    end
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
                
                if sum(absFlag)>0
                    absI = prod( abs( A(r_lag+1:Lt,absFlag)-A(1:Lt-r_lag,absFlag) ) ,2);
                else
                    absI = 1;
                end
                
                I = difI.*midI.*absI; % increments
                
                S(i_s).(var)(r_lag) = mean(I,'omitnan');        % structure function
                
                if iferrors
                    V(i_s).(var)(r_lag) = mean(I.^2,'omitnan');     % variance of increments
                    L(i_s).(var)(r_lag) = int_ls_short(I)*S(i_s).dr;% integral length scale for increments
                    S(i_s).N    (r_lag) = sum(~isnan(I));           % number of valid samples
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
    
    % s3A
    if all(ismember({'uuA','vvA','wwA'},fieldnames(avS)))
        avS(i_l).s3a = 3*( avS(i_l).uuA + avS(i_l).vvA + avS(i_l).wwA );
        avU(i_l).s3a = 3*sqrt( avU(i_l).uuA.^2 + avU(i_l).vvA.^2 + avU(i_l).wwA.^2 );
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
    r = (1:length(avS(i_l).uuu))*avS(i_l).dr;
    
    avS(i_l).W  = 6*cumtrapz(r,avS(i_l).wB .*r.^2) ./ r.^3;
    avS(i_l).Wt = 6*cumtrapz(r,avS(i_l).wBt.*r.^2) ./ r.^3;
    avS(i_l).Wq = 6*cumtrapz(r,avS(i_l).wBq.*r.^2) ./ r.^3;
    avU(i_l).W  = 6*cumtrapz(r,avU(i_l).wB .*r.^2) ./ r.^3;
    avU(i_l).Wt = 6*cumtrapz(r,avU(i_l).wBt.*r.^2) ./ r.^3;
    avU(i_l).Wq = 6*cumtrapz(r,avU(i_l).wBq.*r.^2) ./ r.^3;
    
    avS(i_l).s3lr = avS(i_l).s3l./r;
    avU(i_l).s3lr = avU(i_l).s3l./r;
    
    avS(i_l).Ti = 3*cumtrapz(r,avS(i_l).s2W.*r.^2) ./ r.^3;
    avU(i_l).Ti = 3*cumtrapz(r,avU(i_l).s2W.*r.^2) ./ r.^3;
    
    avS(i_l).Ws3lr = avS(i_l).W - avS(i_l).s3lr;
    avU(i_l).Ws3lr = avU(i_l).W + avU(i_l).s3lr;
    
    avS(i_l).ms3l = -avS(i_l).s3l;
    avU(i_l).ms3l = avU(i_l).s3l;
    
    avS(i_l).Wrs3l = avS(i_l).W.*r - avS(i_l).s3l;
    avU(i_l).Wrs3l = avU(i_l).W.*r + avU(i_l).s3l;
    
    avS(i_l).uuu3r = 3*avS(i_l).uuu./r;
    avS(i_l).vvu3r = 3*avS(i_l).vvu./r;
    avS(i_l).wwu3r = 3*avS(i_l).wwu./r;
    avU(i_l).uuu3r = 3*avU(i_l).uuu./r;
    avU(i_l).vvu3r = 3*avU(i_l).vvu./r;
    avU(i_l).wwu3r = 3*avU(i_l).wwu./r;
    
    
    if all(ismember({'uuA','vvA','wwA'},fieldnames(avS)))
        avS(i_l).s3ar = avS(i_l).s3a./r;
        avU(i_l).s3ar = avU(i_l).s3a./r;

        avS(i_l).Ws3ar = avS(i_l).W - avS(i_l).s3ar;
        avU(i_l).Ws3ar = avU(i_l).W + avU(i_l).s3ar;

        avS(i_l).uuA3r = 3*avS(i_l).uuA./r;
        avS(i_l).vvA3r = 3*avS(i_l).vvA./r;
        avS(i_l).wwA3r = 3*avS(i_l).wwA./r;
        avU(i_l).uuA3r = 3*avU(i_l).uuA./r;
        avU(i_l).vvA3r = 3*avU(i_l).vvA./r;
        avU(i_l).wwA3r = 3*avU(i_l).wwA./r;
    end
    
end



%% Estimate transport

% (1 idea) as residual

% Nlvl = size(avS,1);
% for i_l = 1:Nlvl
%     avS(i_l).Tres = -avS(i_l).s3lr + avS(i_l).W - 4*avS(i_l).edr_s2;
%     avU(i_l).Tres = sqrt( avU(i_l).s3lr.^2 + avU(i_l).W.^2 + 4*avU(i_l).edr_s2.^2 );
% end


% (2 idea) from differences between levels

% level_pairs = {'cloud-base','top-subcloud';
%                'cloud-base-noclouds','top-subcloud';
%                'top-subcloud','mid-subcloud';
%                'mid-subcloud','near-surface'};
% 
% for i_l = 1:size(level_pairs,1)
%     i1 = find([avS(:).level]==level_pairs{i_l,1});
%     i2 = find([avS(:).level]==level_pairs{i_l,2});
%     
%     avS(i1).Tdif = (avS(i1).Ti-avS(i2).Ti) / (avS(i1).alt-avS(i2).alt);
%     avU(i1).Tdif = sqrt(avU(i1).Ti.^2+avU(i2).Ti.^2) / (avS(i1).alt-avS(i2).alt);
% end


% (3 idea) from interpolated smooth function

alt_increment = 25; % m
int_method = 'pchip';


int_levels = {'cloud-base','top-subcloud','mid-subcloud','near-surface'};
[~,ind_l] = ismember(int_levels,[avS(:).level]);

alt_vec   = [avS(ind_l).alt];
alt_query = reshape( [alt_vec-alt_increment; alt_vec+alt_increment] ,1,[]);
% alt_query = reshape( [alt_vec; alt_vec+alt_increment] ,1,[]);

avS_T_int = interp1(alt_vec,cat(1,avS(ind_l).Ti),alt_query,int_method);
if iferrors
    avU_T_int = interp1(alt_vec,cat(1,avU(ind_l).Ti),alt_query,int_method);
else
    avU_T_int = nan(size(avS_T_int));
end

for i_l = 1:numel(ind_l)
    avS(ind_l(i_l)).Tint = ( avS_T_int(2*i_l,:) - avS_T_int(2*i_l-1,:) )/alt_increment/2;
    avU(ind_l(i_l)).Tint = sqrt( avU_T_int(2*i_l,:).^2 + avU_T_int(2*i_l-1,:).^2 )/alt_increment/2;
end


int_levels = {'cloud-base-noclouds','top-subcloud','mid-subcloud','near-surface'};
[~,ind_l] = ismember(int_levels,[avS(:).level]);

alt_vec   = [avS(ind_l).alt];
alt_query = reshape( [alt_vec-alt_increment; alt_vec+alt_increment] ,1,[]);
% alt_query = reshape( [alt_vec; alt_vec+alt_increment] ,1,[]);

avS_T_int = interp1(alt_vec,cat(1,avS(ind_l).Ti),alt_query,int_method);
if iferrors
    avU_T_int = interp1(alt_vec,cat(1,avU(ind_l).Ti),alt_query,int_method);
else
    avU_T_int = nan(size(avS_T_int));
end

i_l = 1;
avS(ind_l(i_l)).Tint = ( avS_T_int(2*i_l,:) - avS_T_int(2*i_l-1,:) )/alt_increment/2;
avU(ind_l(i_l)).Tint = sqrt( avU_T_int(2*i_l,:).^2 + avU_T_int(2*i_l-1,:).^2 )/alt_increment/2;


% Dependent parameters

for i_l = 1:Nlvl
    r = (1:length(avS(i_l).uuu))*avS(i_l).dr;

    avS(i_l).WST = avS(i_l).W - avS(i_l).s3lr - avS(i_l).Tint;
    avU(i_l).WST = avU(i_l).W + avU(i_l).s3lr + avU(i_l).Tint;
    
    avS(i_l).WrST = avS(i_l).W.*r - avS(i_l).s3l - avS(i_l).Tint.*r;
    avU(i_l).WrST = avU(i_l).W.*r + avU(i_l).s3l + avU(i_l).Tint.*r;
end



%% Calculate dissipation rate

fit_range = [20 60];

edr_vars = {'uu','vv','ww','ms3l','Wrs3l','WrST'};
slps = [2/3 2/3 2/3 1 1 1];
cons = [2.0 2.6 2.6 4.0 4.0 4.0];

Nvar = numel(edr_vars);
for i_v = 1:Nvar
    var = edr_vars{i_v};
    [avS,avU] = calc_edr(avS,avU,edr_vars{i_v},'Slope',slps(i_v),'Constant',cons(i_v),...
        'FitRange',fit_range,'Method','direct','FitPoints',6);
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
    avS(i_l).edr_s2_4 = avS(i_l).edr_s2*4;
    avU(i_l).edr_s2_4 = avU(i_l).edr_s2*4;
end



%% Save/load

save(['S_eureca_1km_',orientation,'.mat'],'S','U','N','L','V','avS','avU','avN',...
    'r_max','dr','r_maxlag','r','MOM','levels','plotpath',...
    'MOM2','sfc_vars','valid_seg')


% load(['S_eureca_1km_',orientation,'.mat'])
% 
% addpath(genpath(myprojectpath))
% plotpath = [myprojectpath,filesep,'figures_',orientation];
% if ~isfolder(plotpath)
%     mkdir(plotpath)
% end



%% TABLES

% Segment overview

print_table(MOM2(valid_seg,:),{'length_km','alt'},1,0)


% Dissipation rates

print_vars = {'uu','vv','ww','s2','ms3l','Wrs3l','WrST'};

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


%% Overview of the segments

plot_seg_overview(MOM,setdiff(levels,{'cloud-base-noclouds'},'stable'),false);
print(gcf,[plotpath,filesep,'seg_overview'],'-dpng','-r300')


%% Structure functions

% plots = {
%     {'uuu3r','vvu3r','wwu3r'}, {'uuu','vvu','wwu'}, '$[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
%     {'W','Wt','Wq'}, {'$W$','$W_\theta$','$W_q$'},'$W\,[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
%     {'s3lr'}, {'$S_3^L r^{-1}$'}, '$S_3^Lr^{-1}\,[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
%     {'Ws3lr','edr_s2_4'}, {'$W-S_3^L r^{-1}$','$4\epsilon_2$'},'$[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
%     {'Tres','Tdif'}, {'$T_{res}$','$T_{dif}$'}, '$T\,[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
%     {'uu_c','vv_c','ww_c'}, {'uu','vv','ww'}, '$S_2r^{-2/3}C^{-1}$', [8e-4 1e-2];
%     {'W','m_s3lr','m_Tdif','edr_s2_4'}, {'$W$','$-S_3^L r^{-1}$','$-T_{dif}$','$4\epsilon_2$'}, '$[\mathrm{m^2\,s^{-3}}]$', [1e-6 4e-3];
%     };
% %%

xlim = [4 1000];
ylim = [1e-6 5e-3];
% ylim = [1e-6 1e-1];

Npoints = 40;

Nlvl = size(avS,1);


for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','Wt','Wq'},[7 5 1],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l==1
        legend({'$W$','$W_\theta$','$W_q$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'W_',levels{i_l}],'-dpng','-r300')
end

%%
for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'s3lr'},[1 6],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l==1
        legend({'$S_3 r^{-1}$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s3l_',levels{i_l}],'-dpng','-r300')
end

%%
for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'WST','edr_s2_4'},[4 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l==1
        legend({'$W-S_3 r^{-1}-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt)) 
    print(fig,[plotpath,filesep,'WST_',levels{i_l}],'-dpng','-r300')
end

%%
for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'W','-s3lr','-Tint','edr_s2_4'},[7 1 5 3],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l==1
        legend({'$W$','$-S_3 r^{-1}$','$-T_u$','$4\epsilon_2$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'balance_',levels{i_l}],'-dpng','-r300')
end

%%
for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'uu_c','vv_c','ww_c'},[7 5 1],Npoints,'XLim',xlim,'YLim',[8e-4 1e-2]);
    if i_l==1
        legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    end
    ylabel('$S_2r^{-2/3}C^{-1} $','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s2_',levels{i_l}],'-dpng','-r300')
end

%%
for i_l = 1:Nlvl
    fig = plot_sfc(avS(i_l),avU(i_l),{'uuu3r','vvu3r','wwu3r'},[7 5 1 2 4 6],Npoints,'XLim',xlim,'YLim',ylim);
    if i_l==1
        legend({'uuu','vvu','wwu'},'Interpreter','latex','Location','best')
    end
    ylabel('$3S_3 r^{-1}\,[\mathrm{m^2\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'s3_comp_',levels{i_l}],'-dpng','-r300')
end



% for i_p = 1:size(plots,1)
%     for i_l = 1:Nlvl
%         fig = plot_sfc(avS(i_l),avU(i_l),plots{i_p,1},[7 1 5 3],...
%             'XLim',xlim,'YLim',plots{i_p,4});
%         if numel(plots{i_p,2})>1
%             legend(plots{i_p,2},'Interpreter','latex','Location','best')
%         end
%         ylabel(plots{i_p,3},'Interpreter','latex')
%         title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
%         print(fig,[plotpath,filesep,plots{i_p,1}{1},'_',levels{i_l}],'-dpng','-r300')
%     end
% end



%% Fitting dissipation

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    
    fig = plot_sfc_edr(avS(i_l),[],{'uu','vv','ww'},[7 5 1],Npoints,'XLim',xlim,'YLim',[1e-3 2e-2]);
    if i_l==1
        legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    end
    ylabel('$S_2r^{-2/3}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr2_',levels{i_l}],'-dpng','-r300')
end

%%
for i_l = 1:Nlvl
    
    fig = plot_sfc_edr(avS(i_l),[],{'ms3l','Wrs3l','WrST'},[1 7 4],Npoints,'XLim',xlim,'YLim',[1e-5 1e-2]);
    if i_l==1
        legend({'$-S_3$','$-S_3+Wr$','$-S_3+Wr-T$'},'Interpreter','latex','Location','best')
    end
    ylabel('$[\mathrm{m^3\,s^{-3}}]$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'edr3_',levels{i_l}],'-dpng','-r300')
    
end
