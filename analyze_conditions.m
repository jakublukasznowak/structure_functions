
%% Settings

% Cloud mask

mask_rh_threshold = 98;   % Relative humidity threshold for cloud mask
mask_ext_factor   = 1;    % Mask extension factor to account for cloud shells
mask_nan_replace  = true; % Assume cloud if masks are nan



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

% List of variables from core dataset
core_vars = {'time','ALTITUDE','TAS','HEADING','WIND_FF','WIND_DD','TEMPERATURE','LWC'};

           

%% Prepare paths

addpath(genpath(myprojectpath))

datapath = mydatapath;

plotpath = [myprojectpath,filesep,'figures',filesep,'conditions'];
if ~isfolder(plotpath), mkdir(plotpath), end



%% Load datasets

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

MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';


% Microphysics

disp('Load microphysics data ...')
MCPH = load_atr_pma(MOM,datapath,mcph_vars);


% Core variables

disp('Load core data ...')
CORE = load_atr_core(MOM,datapath,core_vars);


Nseg = size(MOM,1);
fprintf('Number of segments loaded: %d\n',Nseg)



%% Compute auxilliary parameters

disp('Compute thermodynamics and cloud mask ...')


% Thermodynamics (including buoyancy)

TURB = calc_thermo(TURB,MOM);


% Cloud masks

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



%% Compute fluxes

disp('Compute fluxes ...')

Nseg2 = size(MOM,1);

for i_s = 1:Nseg2
    w = TURB(i_s).W_DET;
    if MOM.level(i_s) == "cloud-base-noclouds"
        w(logical(TURB(i_s).OR_mask)) = nan;
    end
    
    MOM.WU(i_s) = mean( w.*TURB(i_s).UX_DET, 'omitnan' );
    MOM.WV(i_s) = mean( w.*TURB(i_s).VY_DET, 'omitnan' );
    MOM.WB(i_s) = mean( w.*TURB(i_s).B, 'omitnan' );
end

MOM.Ustar2 = (MOM.WU.^2 + MOM.WV.^2).^0.5;
MOM.Ustar = (MOM.WU.^2 + MOM.WV.^2).^0.25;
k = 0.4;
MOM.OB = -MOM.Ustar.^3/k./MOM.WB;
MOM.zL = MOM.alt./MOM.OB;



%% PLOTS

lw = 1;

[f,~,co] = fig16x12;
close(f)
co = repmat(co,100,1);


%% Fluxes

h = plot_whisker(MOM,{'WU','WV','Ustar2'},levels,0,...
    'PrimaryLabels',{'$\langle wu \rangle$','$\langle wv \rangle$','$u*^2$'});
hold on
ylabel('$[\mathrm{m^2s^{-2}}]$','Interpreter','latex')
print(h.figure,[plotpath,filesep,'flux_momentum'],'-dpng','-r300')

h = plot_whisker(MOM,{'WB'},levels,0,'PrimaryLabels',{''});
hold on
ylabel('$\langle wB \rangle\,[\mathrm{m^2s^{-2}}]$','Interpreter','latex')
print(h.figure,[plotpath,filesep,'flux_wB'],'-dpng','-r300')

h = plot_whisker(MOM,{'OB'},levels,0,'PrimaryLabels',{''},...
    'Whisker',2,'datalim',[-1000,1000],'extrememode','clip');
hold on
ylabel('$Ob\,[\mathrm{m}]$','Interpreter','latex')
print(h.figure,[plotpath,filesep,'flux_Ob'],'-dpng','-r300')

h = plot_whisker(MOM,{'zL'},levels,0,'PrimaryLabels',{''},...
    'Whisker',2,'datalim',[-10,10],'extrememode','clip');
hold on
ylabel('$z/Ob$','Interpreter','latex')
print(h.figure,[plotpath,filesep,'flux_zL'],'-dpng','-r300')


%% Basic meteo

levels = setdiff(levels,{'cloud-base-noclouds'});

Nlvl = numel(levels);
Nvar = numel(core_vars);
Nseg = size(CORE,1);

for i_v = 2:Nvar
    var = core_vars{i_v};
    
    for i_l = 1:Nlvl
        lvl = levels{i_l};
    
        fig = figure('Units','normalized','Position',[0 0.5 1 0.4]);
        grid on, hold on
        
        i_c = 1; legstr = {};
        
        for i_s = 1:Nseg
            time0 = CORE(i_s).time(1);
            time  = seconds(CORE(i_s).time -time0);
            tas = MOM.MEAN_TAS(i_s)/1000; % [km/s]
            xlim([0 tas*time(end)]);
            
            if CORE(i_s).level==lvl
                plot(tas*time,CORE(i_s).(var),'Color',co(i_c,:),'LineWidth',lw);
                i_c = i_c + 1;
                legstr = vertcat(legstr,{CORE(i_s).flight+" "+CORE(i_s).name});
            end
        end
        
        xlabel('Distance [km]')
        ylabel(replace(var,'_','\_'))
        legend(legstr,'Location','eastoutside')
        title(lvl)
        print(fig,[plotpath,filesep,lower(var),'_',lvl],'-dpng','-r300')

    end
    
 end

  
%% Turbulent fluctuations

levels = setdiff(levels,{'cloud-base-noclouds'});

Nvar = numel(turb_vars);
Nlvl = numel(levels);
Nseg = size(TURB,1);

for i_v = 2:Nvar
    var = turb_vars{i_v};
    
    for i_l = 1:Nlvl
        lvl = levels{i_l};
    
        fig = figure('Units','normalized','Position',[0 0 1 0.4]);
        grid on, hold on
        
        i_c = 1; legstr = {};
        
        for i_s = 1:Nseg
            time0 = TURB(i_s).time(1);
            time  = seconds(TURB(i_s).time -time0);
            tas = MOM.MEAN_TAS(i_s)/1000; % [km/s]
            xlim([0 tas*time(end)]);
            
            if TURB(i_s).level==lvl
                plot(tas*time,TURB(i_s).(var),'Color',co(i_c,:),'LineWidth',lw);
                i_c = i_c + 1;
                legstr = vertcat(legstr,{TURB(i_s).flight+" "+TURB(i_s).name});
            end
        end
        
        xlabel('Distance [km]')
        ylabel(replace(var,'_','\_'))
        legend(legstr,'Location','eastoutside')
        title(lvl)
        print(fig,[plotpath,filesep,lower(var),'_',lvl],'-dpng','-r300')

    end
    
 end

