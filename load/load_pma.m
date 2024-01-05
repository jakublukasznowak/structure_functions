
function [PMA,pma_info] = load_pma(MOM,datapath,pma_vars)

epoch = datetime('2020-01-01 00:00:00.000');

if ~ismember('time',pma_vars)
    pma_vars = cat(2,{'time'},pma_vars);
end


% Read nc files

d = dir([datapath,filesep,'CloudComposite',filesep,'*.nc']);

Nf = numel(d);
pma = cell(Nf,1);
pma_info = cell(Nf,1);

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
    PMA(i_s).name   = MOM.name(i_s);
    PMA(i_s).level  = MOM.level(i_s);
    PMA(i_s).type   = MOM.type(i_s);
    
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

end
