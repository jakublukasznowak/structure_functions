
% DATAPATH is the path where you downloaded two datasets:
%
%
% DATAPATH/TURBLENCE
%
% Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz. Aeris.
% https://doi.org/10.25326/128
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/
% 
% In this code only 'longlegs' L3 v1.9 is used.
%
% 
% DATAPATH/CloudComposite
%
% Coutris, P. (2021). SAFIRE ATR42: PMA/Cloud composite dataset. Aeris.
% https://doi.org/10.25326/237
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/PMA/PROCESSED/CloudComposite/
% 
% In this code only v1 is used.


addpath([pwd,filesep,'load'])
addpath([pwd,filesep,'utils'])
addpath([pwd,filesep,'utils',filesep,'yaml_Koch'])


% List of flights
flight_ids = num2cell(num2str((9:19)','RF%02d'),2); % RF09 - RF19

% List of variables from turbulent moments dataset
mom_vars = {'alt','time_start','time_end',...
    'MEAN_WSPD','MEAN_WDIR','MEAN_TAS','MEAN_THDG'};

% List of variables from turbulent fluctuations dataset
turb_vars = {'time','W_DET','UX_DET','VY_DET','T_DET','MR_DET'};

% List of variables from cloud composite dataset
pma_vars = {'time','LWC','NT','CLOUD_mask'};




SEG = load_seg(mydatapath,'v1.9','longlegs');
SEG = SEG(ismember(SEG.flight,flight_ids),:);

[MOM,mom_info] = load_mom(mydatapath,'L3','v1.9','longlegs',mom_vars);
MOM = join(SEG,MOM,'Keys',{'start','xEnd'});

TURB = load_turb(MOM,mydatapath,'L3','v1.9',turb_vars);

PMA = load_pma(MOM,mydatapath,pma_vars);



%% Checkout levels

levels  = {'cloud-base','top-subcloud','mid-subcloud','near-surface'};


fig = figure('Units','normalized','Position',[0 0 0.6 0.3]);
hold on, grid on
co = get(gca,'ColorOrder');

Nlvl = numel(levels);
for i_l = 1:Nlvl
    ind_l = (MOM.level==levels{i_l});
    plot(MOM{ind_l,"start"}, MOM{ind_l,"alt"},...
        'Marker','o','MarkerSize',10,'LineStyle','none',...
        'Color',co(i_l,:),'MarkerFaceColor',co(i_l,:))
    text(MOM{ind_l,"start"}, MOM{ind_l,"alt"}, MOM{ind_l,"name"})
end

legend(levels)
ylabel('Altitude [m]')

