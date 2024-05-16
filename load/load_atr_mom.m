
function [MOM,mom_info] = load_atr_mom(datapath,datalevel,dataversion,segtype,vars)

if nargin<5
    vars = {};
end
if ~isempty(vars) && ~all(ismember({'time_start','time_end'},vars))
    vars = union({'time_start','time_end'},vars,'stable');
end


d = dir([datapath,filesep,'TURBULENCE',filesep,'TURB_MOMENTS',filesep,...
    datalevel,filesep,dataversion,filesep,segtype,filesep,'*.nc']);
Nf = numel(d);


MOM = cell(Nf,1); mom_info = cell(Nf,1);
for i_f = 1:Nf
%     fprintf('Load %s\n',d(i_f).name)
    [MOM{i_f},mom_info{i_f}] = load_nc([d(i_f).folder,filesep,d(i_f).name],vars);
    MOM{i_f} = struct2table(MOM{i_f});
end

MOM = cat(1,MOM{:});
mom_info = cat(1,mom_info{:});


epoch = datetime('2020-01-01 00:00:00.000');
MOM.start = datetime(MOM.time_start,'ConvertFrom','epochtime',...
    'Epoch',epoch,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SS');
MOM.end   = datetime(MOM.time_end,  'ConvertFrom','epochtime',...
    'Epoch',epoch,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SS');
MOM(:,{'time_end','time_start'}) = [];

MOM.alt = double(MOM.alt);
MOM.length = seconds(MOM.end-MOM.start).*MOM.MEAN_TAS;

if all(ismember({'MEAN_THDG','MEAN_WDIR'},MOM.Properties.VariableNames))
    [MOM.dir2,MOM.dir4] = dev_angle(MOM.MEAN_THDG,MOM.MEAN_WDIR);
end


end