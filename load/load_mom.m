
function [MOM,mom_info] = load_mom(datapath,datalevel,dataversion,segtype,mom_vars)

if ~ismember('time',mom_vars)
    mom_vars = cat(2,{'time'},mom_vars);
end


d = dir([datapath,filesep,'TURBULENCE',filesep,'TURB_MOMENTS',filesep,...
    datalevel,filesep,dataversion,filesep,segtype,filesep,'*.nc']);
Nf = numel(d);


MOM = cell(Nf,1); mom_info = cell(Nf,1);
for i_f = 1:Nf
    fprintf('Load %s\n',d(i_f).name)
    [MOM{i_f},mom_info{i_f}] = load_nc([d(i_f).folder,filesep,d(i_f).name],mom_vars);
    MOM{i_f} = struct2table(MOM{i_f});
end

MOM = cat(1,MOM{:});
mom_info = cat(1,mom_info{:});


epoch = datetime('2020-01-01 00:00:00.000');
MOM.start = datetime(MOM.time_start,'ConvertFrom','epochtime',...
    'Epoch',epoch,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SS');
MOM.xEnd   = datetime(MOM.time_end,  'ConvertFrom','epochtime',...
    'Epoch',epoch,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SS');


MOM.alt = double(MOM.alt);

end