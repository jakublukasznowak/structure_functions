
function SEG = load_atr_seg (datapath,dataversion,segtype)

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


d = dir([datapath,filesep,'TURBULENCE',filesep,'YAML',filesep,...
    dataversion,filesep,'yaml_',segtype,filesep,'*.yaml']);
Nf = numel(d);

SEG = cell(Nf,1);
for i_f = 1:Nf
    fprintf('Load %s\n',d(i_f).name)
    yml = yaml.loadFile([d(i_f).folder,filesep,d(i_f).name],'ConvertToArray',true);
    for i_s = 1:numel(yml.legs)
        yml.legs(i_s).flight = yml.flight_id;
        yml.legs(i_s).level = level_dict(strcmp(level_dict(:,1),yml.legs(i_s).kind),2);
        yml.legs(i_s).type = string(segtype);
    end
    SEG{i_f} = yml.legs';
end

SEG = cat(1,SEG{:});
SEG = struct2table(SEG);

SEG.end = SEG.xEnd;
SEG.xEnd = [];

frontfields = {'flight','name','level','type','start','end'};
SEG = movevars(SEG,frontfields,'Before',1);

end
