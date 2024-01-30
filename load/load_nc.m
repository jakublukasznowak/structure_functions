
function [data,info] = load_nc(file,vars,varsnew)

monit='UNSUPPORTED DATATYPE';

if nargin<3
    varsnew = {};
end
if nargin<2
    vars = {};
end


info = ncinfo(file);

if isempty(vars)
    vars = {info.Variables.Name};
end
if isempty(varsnew)
    varsnew = vars;
end


for i = 1:numel(vars)
    
    [if_var,ind_var] = ismember(vars{i},{info.Variables.Name});
    
    if if_var
        if strcmp(info.Variables(ind_var).Datatype,monit)
            fprintf('Warning in LOAD_NC: Variable %s is of %s.\n',vars{i},monit)
        else
            ind = find(~isstrprop(varsnew{i},'digit'),1,'first');
            data.(varsnew{i}(ind:end)) = ncread(file,vars{i});
            if ind>1
                info.Variables(ind_var).Name = info.Variables(ind_var).Name(ind:end);
            end
        end
    else
        fprintf('Warning in LOAD_NC: Variable %s not found in the file.\n',vars{i})
    end
    
end
    
end