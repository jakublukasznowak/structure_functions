
addpath('load','calculate','cloud_mask','plots','utils','utils/yaml_Koch')



%% Load datasets

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

SEG = load_seg(mydatapath,'v1.9','longlegs');
SEG = SEG(ismember(SEG.flight,flight_ids) & ismember(SEG.level,levels),:);

[MOM,mom_info] = load_mom(mydatapath,'L3','v1.9','longlegs',mom_vars);
MOM = join(SEG,MOM,'Keys',{'start','xEnd'});

TURB = load_turb(MOM,mydatapath,'L3','v1.9',turb_vars);

PMA = load_pma(MOM,mydatapath,pma_vars);


% Calculate auxiliary parameters

[MOM.dir2,MOM.dir4] = dev_angle(MOM.MEAN_THDG,MOM.MEAN_WDIR);
MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';

TURB = calculate_thermo(TURB,MOM); % thermodynamic parameters


% Cloud masks
% (1) from PMA dataset based on LWC and interpolated to 25 Hz
% (2) based on RH>98%
% both extended in front and behind each cloud by its 1 width
% (3) any of the two extended masks

[TURB,MOM] = cloud_mask(TURB,MOM,PMA,98,1);


% Plot overview of the segments

plot_seg_overview(MOM);


















