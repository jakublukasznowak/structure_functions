
function SEG = load_seg (datapath,dataversion,segtype)


d = dir([datapath,filesep,'TURBULENCE',filesep,'YAML',filesep,...
    dataversion,filesep,'yaml_',segtype,filesep,'*.yaml']);
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
        yml.legs(i_s).type = string(segtype);
    end
    SEG{i_f} = yml.legs';
end

SEG = cat(1,SEG{:});
SEG = struct2table(orderfields(SEG,{'flight','name','level','type','kind','start','xEnd'}));

end
