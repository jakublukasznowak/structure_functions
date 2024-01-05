
% TODO
% - correct number of points and valid fractions


projectpath = pwd;

plotpath = [projectpath,filesep,'figures'];

datapath = mydatapath;


% DATAPATH is the path where you downloaded two datasets:
%
% DATAPATH/TURBLENCE
%
% Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz. Aeris.
% https://doi.org/10.25326/128
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/
% 
% In this code 'longlegs' L3 v1.9 is used.
% 
% DATAPATH/CloudComposite
%
% Coutris, P. (2021). SAFIRE ATR42: PMA/Cloud composite dataset. Aeris.
% https://doi.org/10.25326/237
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/PMA/PROCESSED/CloudComposite/
% 
% In this code v1 is used.


% Prepare paths

addpath(genpath(projectpath))
if ~isfolder(plotpath)
    mkdir(plotpath)
end



%% Load datasets

% List of flights

flight_ids = num2cell(num2str((9:19)','RF%02d'),2); % RF09 - RF19

% List of levels

levels  = {'cloud-base','top-subcloud','mid-subcloud','near-surface'};

% List of variables from turbulent moments dataset

mom_vars = {'alt','time_start','time_end',...
    'MEAN_P','MEAN_THETA','MEAN_MR',...
    'MEAN_WSPD','MEAN_WDIR','MEAN_TAS','MEAN_THDG'};

% List of variables from turbulent fluctuations dataset

turb_vars = {'time','W_DET','UX_DET','VY_DET','T_DET','MR_DET'};

% List of variables from cloud composite dataset

pma_vars = {'time','LWC','NT','CLOUD_mask'};


% Read data files

SEG = load_seg(datapath,'v1.9','longlegs');
SEG = SEG(ismember(SEG.flight,flight_ids) & ismember(SEG.level,levels),:);

[MOM,mom_info] = load_mom(datapath,'L3','v1.9','longlegs',mom_vars);
MOM = join(SEG,MOM,'Keys',{'start','xEnd'});

TURB = load_turb(MOM,datapath,'L3','v1.9',turb_vars);

PMA = load_pma(MOM,datapath,pma_vars);


% Calculate auxiliary parameters

[MOM.dir2,MOM.dir4] = dev_angle(MOM.MEAN_THDG,MOM.MEAN_WDIR);
MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';
MOM.length = (MOM.time_end-MOM.time_start).*MOM.MEAN_TAS;

TURB = calc_thermo(TURB,MOM); % thermodynamic parameters


% Cloud masks
% (1) from PMA dataset based on LWC and interpolated to 25 Hz
% (2) based on RH>98%
% both extended in front and behind each cloud by its 1 width
% (3) any of the two extended masks

[TURB,MOM] = cloud_mask(TURB,MOM,PMA,98,1);

MOM.valid = 1 - MOM.OR_cloud_fraction;
MOM.valid(MOM.level~="cloud-base") = 1;


% Plot overview of the segments

plot_seg_overview(MOM);



%% Calculate structure functions

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
            'uuW',{'VY_DET','VY_DET','W_MID'};
            'vvW',{'UX_DET','UX_DET','W_MID'};
            'wwW',{'W_DET','W_DET','W_MID'}
            };

% Assume fixed distance between points
MOM.dr(:) = 4;

% Consider max lag of 400 m
r_maxlag = 400/4;


% Calculate velocity in between data points
Nseg = size(TURB,1);
for i_s = 1:Nseg
    Lt = length(TURB(i_s).time);
    TURB(i_s).W_MID = interp1(TURB(i_s).W_DET,0.5:0.5:Lt)';
    TURB(i_s).OR_mask_ext_mid = (interp1(TURB(i_s).W_DET,0.5:0.5:Lt)'>0);
end
    


S = table2struct(MOM(:,{'flight','name','level','alt','dr','length','valid'}));
V = S; L = S; N = S; U = S;

Nseg = size(S,1);
Nvar = size(sfc_vars,1);

for i_v = 1:Nvar
    
    var = sfc_vars{i_v,1};
    fprintf('S_%s = < ',var)
    fprintf('d %s ',sfc_vars{i_v,2}{:})
    fprintf('>\n')
    
    for i_s = 1:Nseg
        
        Lt = length(TURB(i_s).time);
    
        A = cellfun(@(v) TURB(i_s).(v),sfc_vars{i_v,2},'UniformOutput',false);
        La = cellfun(@length,A);
        
        difA = cell2mat(A(La==Lt));
        midA = cell2mat(A(La>Lt));
        
        if S(i_s).level == "cloud-base"
            difA(TURB(i_s).OR_mask_ext,:) = nan; % mask cloudy points
            if ~isempty(midA)
                midA(TURB(i_s).OR_mask_ext_mid,:) = nan; % mask cloudy points
            end
        end
        
        if ~any(cellfun(@(x) all(isnan(x),1),A))
            
            S(i_s).(var) = nan(1,r_maxlag);
            V(i_s).(var) = nan(1,r_maxlag);
            L(i_s).(var) = nan(1,r_maxlag);
            S(i_s).N     = nan(1,r_maxlag);

            for r_lag = 1:r_maxlag
                
                difI = prod( difA(r_lag+1:Lt,:)-difA(1:Lt-r_lag,:) ,2);
                
                if ~isempty(midA)
                    midI = prod( midA(r_lag+1:2:2*Lt-r_lag,:) ,2);    
                else
                    midI = 1;
                end
                
                I = difI.*midI; % increments
                
                S(i_s).(var)(r_lag) = mean(I,'omitnan');        % structure function
                V(i_s).(var)(r_lag) = mean(I.^2,'omitnan');     % variance of increments
                L(i_s).(var)(r_lag) = int_ls_short(I)*S(i_s).dr;% integral length scale for increments
                S(i_s).N    (r_lag) = sum(~isnan(I));           % number of increment samples
            end
        
            N(i_s).(var) = S(i_s).valid*S(i_s).length./L(i_s).(var); % number of independent samples
            % N(i_s).(var) = S(i_s).N*dr./L(i_s).(var); % alternative formula
            U(i_s).(var) = sqrt( (V(i_s).(var) - S(i_s).(var).^2) ./ N(i_s).(var) ); % uncertainty
              
        end

    end
    
end



%% Average structure functions at levels

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
valid_hum = (MOM.RH_nan_fraction==0);


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
        avN(i_l,1).(var) = sum( repmat([S(ind_lv).valid]'.*[S(ind_lv).length]',1,r_maxlag) ...
            ./ cat(1,L(ind_lv).(var)) );
        avU(i_l,1).(var) = sqrt( (avV - avS(i_l).(var).^2) ./ avN(i_l).(var) );
    end
    
    % s3L
    if all(ismember({'uuu','vvu','wwu'},fieldnames(avS)))
        avS(i_l).s3l = 3*( avS(i_l).uuu + avS(i_l).vvu + avS(i_l).wwu );
        avU(i_l).s3l = 3*sqrt( avU(i_l).uuu.^2 + avU(i_l).vvu.^2 + avU(i_l).wwu.^2 );
    end
    
    % wBt + wBq
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
    
    avS(i_l).s3lr = avS(i_l).s3l./r;
    avS(i_l).uuu3r = 3*avS(i_l).uuu./r;
    avS(i_l).vvu3r = 3*avS(i_l).vvu./r;
    avS(i_l).wwu3r = 3*avS(i_l).wwu./r;
    
    avU(i_l).s3lr = avU(i_l).s3l./r;
    avU(i_l).uuu3r = 3*avU(i_l).uuu./r;
    avU(i_l).vvu3r = 3*avU(i_l).vvu./r;
    avU(i_l).wwu3r = 3*avU(i_l).wwu./r; 

    avS(i_l).W  = 6*cumtrapz(r,avS(i_l).wB .*r.^2) ./ r.^3;
    avS(i_l).Wt = 6*cumtrapz(r,avS(i_l).wBt.*r.^2) ./ r.^3;
    avS(i_l).Wq = 6*cumtrapz(r,avS(i_l).wBq.*r.^2) ./ r.^3;
    avS(i_l).T  = 3*cumtrapz(r,avS(i_l).s2W.*r.^2) ./ r.^3;
    
    avU(i_l).W  = 6*cumtrapz(r,avU(i_l).wB .*r.^2) ./ r.^3;
    avU(i_l).Wt = 6*cumtrapz(r,avU(i_l).wBt.*r.^2) ./ r.^3;
    avU(i_l).Wq = 6*cumtrapz(r,avU(i_l).wBq.*r.^2) ./ r.^3;
    avU(i_l).T  = 3*cumtrapz(r,avU(i_l).s2W.*r.^2) ./ r.^3;
    
    avS(i_l).s3lrW = avS(i_l).s3lr - avS(i_l).W;
    avU(i_l).s3lrW = avU(i_l).s3lr + avU(i_l).W;
    avS(i_l).Wrs3l = -avS(i_l).s3l + avS(i_l).W.*r;
    avU(i_l).Wrs3l =  avU(i_l).s3l + avU(i_l).W.*r;
end



%% Plot structure functions

plots = { {'uuu3r','vvu3r','wwu3r'}, {'uuu','vvu','wwu'};
          {'Wt','Wq','W'}, {'$W_\theta$','$W_q$','$W$'};
          {'s3lr','W','s3lrW'}, {'$S_3^L r^{-1}$','$W$','$S_3^L r^{-1}-W$'} };


Nlvl = size(avS,1);
      
for i_p = 1:size(plots,1)
    for i_l = 1:Nlvl
        fig = plot_sfc(avS(i_l),avU(i_l),plots{i_p,1},[5 6 7],...
            'XLim',[4 400],'YLim',[1e-6 2e-3]);
        legend(plots{i_p,2},'Interpreter','latex','Location','best')
        title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
        print(fig,[plotpath,filesep,'S_',levels{i_l},'_',plots{i_p,1}{1}],'-dpng','-r300')
    end
end

fig = plot_sfc(avS,avU,{'T'},[1 2 3 4],'XLim',[4 400]);
legend(levels,'Interpreter','latex','Location','best')
title('Integrated transport')
print(fig,[plotpath,filesep,'S_transport'],'-dpng','-r300')



%% Calculate dissipation rate

edr_vars = {'uu','vv','ww','Wrs3l'};
slps = [2/3 2/3 2/3 1];
cons = [2.0 2.6 2.6 4.0];

fit_range = [10 60];
method = "direct";
fit_points = 6;


Nlvl = size(avS,1);
Nvar = numel(edr_vars);

for i_l = 1:Nlvl
    for i_v = 1:Nvar
        var = edr_vars{i_v};
        avS(i_l).(['fit_',var]) = fit_range;
        avS(i_l).(['slp_',var]) = slps(i_v);
        avS(i_l).(['con_',var]) = cons(i_v);
    end
end

[avS,avU] = calc_edr(avS,avU,edr_vars,'Method',method,'FittingPoints',fit_points);

Tedr = struct2table(rmfield(avS,setdiff(fieldnames(avS),...
    {'level','edr_ww','edr_uu','edr_vv','edr_s2','edr_Wrs3l',...
    'slp_ww_free','slp_uu_free','slp_vv_free','slp_Wrs3l_free'})));



%% Plot fitted scaling laws

Nlvl = size(avS,1);

for i_l = 1:Nlvl
    
    fig = plot_sfc_edr(avS(i_l),[],{'uu','vv','ww'},[5 6 7],'XLim',[4 400],'YLim',[1e-3 2e-2]);
    legend({'uu','vv','ww'},'Interpreter','latex','Location','best')
    ylabel('$S_2r^{-2/3}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'E2_',levels{i_l}],'-dpng','-r300')
    
    fig = plot_sfc_edr(avS(i_l),avU(i_l),{'Wrs3l'},7,'XLim',[4 400],'YLim',[1e-6 2e-3]);
    ylabel('$W-S_3^L r^{-1}$','Interpreter','latex')
    title(sprintf('%s ~%.0f m',levels{i_l},avS(i_l).alt))
    print(fig,[plotpath,filesep,'E3_',levels{i_l}],'-dpng','-r300')
    
end
