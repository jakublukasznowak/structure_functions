
function [CORE,core_info] = load_atr_core(MOM,datapath,core_vars)

if ~ismember('time',core_vars)
    core_vars = cat(2,{'time'},core_vars);
end


% Read nc files

d = dir([datapath,filesep,'CORE',filesep,'*.nc']);

Nf = numel(d);
core = cell(Nf,1);
core_info = cell(Nf,1);

for i_f = 1:Nf
%     fprintf('Load %s\n',d(i_f).name)
    
    [core{i_f},core_info{i_f}] = load_nc([d(i_f).folder,filesep,d(i_f).name],core_vars);
    
    core{i_f}.flight = string(core_info{i_f}.Attributes( ...
        strcmp({core_info{i_f}.Attributes(:).Name},'flight_id')).Value);

    epoch = datetime( core_info{i_f}.Attributes( ...
        strcmp({core_info{i_f}.Attributes(:).Name},'flight_date')).Value );
    
    core{i_f}.fsamp = 1/median(diff(core{i_f}.time));
    core{i_f}.time = datetime(core{i_f}.time,'ConvertFrom','epochtime',...
        'Epoch',epoch,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SS');
    core{i_f}.start = core{i_f}.time(1);
end

core = cat(1,core{:});
core_info = cat(1,core_info{:});
    

% Cut into segments

Nseg = size(MOM,1);
CORE = struct([]);

for i_s = 1:Nseg
    
    CORE(i_s).flight = MOM.flight(i_s);
    CORE(i_s).name   = MOM.name(i_s);
    CORE(i_s).level  = MOM.level(i_s);
    CORE(i_s).type   = MOM.type(i_s);
    
    indF = find([core(:).start]<=MOM.start(i_s),1,'last');
    ind1 = find(core(indF).time>=MOM.start(i_s),1,'first');
    ind2 = find(core(indF).time<=MOM.end(i_s),1,'last');
    
    core_info(indF).Variables(strcmp({core_info(indF).Variables(:).Name},'time_bnds')) = [];
    
    for i_v = 1:numel(core_info(indF).Variables)
        
        if isfield(core(indF),core_info(indF).Variables(i_v).Name)
            
            if ismember('time',{core_info(indF).Variables(i_v).Dimensions(:).Name})
                CORE(i_s).(core_info(indF).Variables(i_v).Name) = ...
                    core(indF).(core_info(indF).Variables(i_v).Name)(ind1:ind2,:);
            else
                CORE(i_s).(core_info(indF).Variables(i_v).Name) = ...
                    core(indF).(core_info(indF).Variables(i_v).Name)';
            end
            
        end
        
    end
    
    CORE(i_s).fsamp = core(indF).fsamp; 
    
end


frontfields = {'flight','name','level','fsamp','time'};
CORE = orderfields(CORE',cat(2,frontfields,setdiff(fieldnames(CORE)',frontfields,'stable')));

end
