
function [TURB,turb_info] = load_atr_turb(MOM,datapath,datalevel,dataversion,vars,varsnew)

epoch = datetime('2020-01-01 00:00:00.000');

if nargin<6
    varsnew = {};
end
if nargin<5
    vars = {};
end
if isempty(varsnew)
    varsnew = vars;
end
if ~isempty(vars) && ~ismember('time',vars)
    vars    = union({'time'},vars,'stable');
    varsnew = union({'time'},varsnew,'stable');
end


Nseg = size(MOM,1);
TURB = cell(Nseg,1);
turb_info = cell(Nseg,1);

for i_s = 1:Nseg
    
    folder = strcat(datapath,filesep,'TURBULENCE',filesep,'TURB_FLUCTUATIONS',filesep,...
        datalevel,filesep,dataversion,filesep,MOM.type(i_s),filesep,MOM.flight(i_s));
    file = strcat('EUREC4A_ATR_turbulent_fluctuations_',datestr(MOM.start(i_s),'yyyymmdd'),...
        '_',MOM.flight(i_s),'_',MOM.name(i_s),'_',datalevel,'_',dataversion,'.nc');
%     fprintf('Load %s\n',file)
    
    [TURB{i_s},turb_info{i_s}] = load_nc(strcat(folder,filesep,file),vars,varsnew);
    
    TURB{i_s}.flight = MOM.flight(i_s);
    TURB{i_s}.name   = MOM.name(i_s);
    TURB{i_s}.level  = MOM.level(i_s);
    TURB{i_s}.type   = MOM.type(i_s);
    
    TURB{i_s}.fsamp = 1/median(diff(TURB{i_s}.time));
    TURB{i_s}.time = datetime(TURB{i_s}.time,'ConvertFrom','epochtime',...
        'Epoch',epoch,'Format','yyyy-MM-dd HH:mm:ss.SS','TimeZone','UTC');
    
end

TURB = cat(1,TURB{:});
turb_info = cat(1,turb_info{:});

frontfields = {'flight','name','level','type','fsamp','time'};
TURB = orderfields(TURB,cat(2,frontfields,setdiff(fieldnames(TURB)',frontfields,'stable')));



end