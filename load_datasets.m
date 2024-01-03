
% DATAPATH is the path where you need to download two datasets:
%
%
% DATAPATH/TURBLENCE
%
% Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz. [Dataset].
% Aeris. https://doi.org/10.25326/128
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/
% 
% In this code only 'longlegs' L3 v1.9 is used.
%
% 
% DATAPATH/CloudComposite
%
% Coutris, P. (2021). SAFIRE ATR42: PMA/Cloud composite dataset. [Dataset].
% Aeris. https://doi.org/10.25326/237
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/PMA/PROCESSED/CloudComposite/
% 
% In this code only v1 is used.


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


% Epoch time
epoch = datetime('2020-01-01 00:00:00.000');



%% Segment information

d = dir([datapath,filesep,'TURBULENCE',filesep,'YAML',filesep,...
    'v1.9',filesep,'yaml_longlegs',filesep,'*.yaml']);
Nf = numel(d);

% Level labels
level_dict = [
    "60m",                "near-surface";
    "Close to surface",   "near-surface";
    "Mid Sub Cld Layer",  "mid-subcloud";
    "Mid Sub Cld layer",  "mid-subcloud";
    "Below Cld base",     "top-subcloud";
    "Just above Cld base","cloud-base";
    "Strati Layer",       "cloud-top";
    "Stratiform layer",   "cloud-top"];

SEG = cell(Nf,1);
for i_f = 1:numel(d)
    fprintf('Load %s\n',d(i_f).name)
    yml = yaml.loadFile([d(i_f).folder,filesep,d(i_f).name],'ConvertToArray',true);
    for i_s = 1:numel(yml.legs)
        yml.legs(i_s).flight = yml.flight_id;
        yml.legs(i_s).level = level_dict(strcmp(level_dict(:,1),yml.legs(i_s).kind),2);
    end
    SEG{i_f} = yml.legs';
end
SEG = cat(1,SEG{:});
SEG = struct2table(orderfields(SEG,{'flight','name','level','kind','start','xEnd'}));

SEG = SEG(ismember(SEG.flight,flight_ids),:);



%% Turbulent moments

d = dir([datapath,filesep,'TURBULENCE',filesep,'TURB_MOMENTS',filesep,...
    'L3',filesep,'v1.9',filesep,'longlegs',filesep,'*.nc']);
Nf = numel(d);

MOM = cell(Nf,1); mom_info = cell(Nf,1);
for i_f = 1:Nf
    fprintf('Load %s\n',d(i_f).name)
    [MOM{i_f},mom_info{i_f}] = load_nc([d(i_f).folder,filesep,d(i_f).name],mom_vars);
    MOM{i_f} = struct2table(MOM{i_f});
end
MOM = cat(1,MOM{:});
mom_info = cat(1,mom_info{:});

MOM.start = datetime(MOM.time_start,'ConvertFrom','epochtime',...
    'Epoch',epoch,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SS');
MOM.xEnd   = datetime(MOM.time_end,  'ConvertFrom','epochtime',...
    'Epoch',epoch,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SS');

MOM = join(SEG,MOM,'Keys',{'start','xEnd'});
MOM.alt = double(MOM.alt);



%% Turbulent fluctuations

Nseg = size(MOM,1);
TURB = cell(Nseg,1); turb_info = cell(Nseg,1);
for i_s = 1:Nseg
    folder = strcat(datapath,filesep,'TURBULENCE',filesep,'TURB_FLUCTUATIONS',filesep,...
        'L3',filesep,'v1.9',filesep,'longlegs',filesep,MOM.flight(i_s));
    file = strcat('EUREC4A_ATR_turbulent_fluctuations_',datestr(MOM.start(i_s),'yyyymmdd'),...
        '_',MOM.flight(i_s),'_',MOM.name(i_s),'_L3_v1.9.nc');
    fprintf('Load %s\n',file)
    
    [TURB{i_s},turb_info{i_s}] = load_nc(strcat(folder,filesep,file),turb_vars);
    
    TURB{i_s}.flight = MOM.flight(i_s);
    TURB{i_s}.name = MOM.name(i_s);
    TURB{i_s}.level = MOM.level(i_s);
    
    TURB{i_s}.fsamp = 1/median(diff(TURB{i_s}.time));
    TURB{i_s}.time = datetime(TURB{i_s}.time,'ConvertFrom','epochtime',...
        'Epoch',epoch,'Format','yyyy-MM-dd HH:mm:ss.SS','TimeZone','UTC');
end
TURB = cat(1,TURB{:});
turb_info = cat(1,turb_info{:});

frontfields = {'flight','name','level','fsamp','time'};
TURB = orderfields(TURB,cat(2,frontfields,setdiff(fieldnames(TURB)',frontfields,'stable')));



%% Cloud composite

% Read nc files

d = dir([datapath,filesep,'CloudComposite',filesep,'*.nc']);
Nf = numel(d);

pma = cell(Nf,1); pma_info = cell(Nf,1);
for i_f = 1:Nf
    fprintf('Load %s\n',d(i_f).name)
    [pma{i_f},pma_info{i_f}] = load_nc([d(i_f).folder,filesep,d(i_f).name],pma_vars);
    pma{i_f}.flight = string(pma_info{i_f}.Attributes( ...
        strcmp({pma_info{i_f}.Attributes(:).Name},'flight_id')).Value);
    pma{i_f}.fsamp = 1/median(diff(pma{i_f}.time));
    pma{i_f}.time = datetime(pma{i_f}.time,'ConvertFrom','epochtime',...
        'Epoch',epoch,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SS');
end
pma = cat(1,pma{:});
pma_info = cat(1,pma_info{:});


% Cut into segments

Nseg = size(MOM,1);
PMA = struct([]);
for i_s = 1:Nseg
    PMA(i_s).flight = MOM.flight(i_s);
    PMA(i_s).name = MOM.name(i_s);
    PMA(i_s).level = MOM.level(i_s);
    
    indF = find([pma(:).flight]==MOM.flight(i_s));
    ind1 = find(pma(indF).time>=MOM.start(i_s),1,'first');
    ind2 = find(pma(indF).time<=MOM.xEnd(i_s),1,'last');
    
    for i_v = 1:numel(pma_info(indF).Variables)
        if isfield(pma(indF),pma_info(indF).Variables(i_v).Name)
            if ismember('time',{pma_info(indF).Variables(i_v).Dimensions(:).Name})
                PMA(i_s).(pma_info(indF).Variables(i_v).Name) = ...
                    pma(indF).(pma_info(indF).Variables(i_v).Name)(ind1:ind2,:);
            else
                PMA(i_s).(pma_info(indF).Variables(i_v).Name) = ...
                    pma(indF).(pma_info(indF).Variables(i_v).Name)';
            end
        end
    end
    PMA(i_s).fsamp = pma(indF).fsamp; 
end

frontfields = {'flight','name','level','fsamp','time'};
PMA = orderfields(PMA',cat(2,frontfields,setdiff(fieldnames(PMA)',frontfields,'stable')));



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

