function tab=extractInstrument(tab, vars, varNames, event, subjs)

if ismember('src_subject_id', tab.Properties.VariableNames)
    tab.src_subject_id = erase(tab.src_subject_id,'_');
end 
if ~isempty(event)
    tab.eventname = categorical(tab.eventname);
    tab = tab(tab.eventname == event,:);
end 
if ~strcmp(vars, 'all')
    tab = tab(:, vars);
end 
if ~isempty(varNames)
    tab = renamevars(tab,vars,varNames);
end

if ~strcmp(subjs, 'all')
    tab = tab(ismember(tab.subjectkey, subjs),:);
end 

for iCol = 1:width(tab)
    currentColumn = tab.(iCol);
    if iscell(currentColumn)
        tab.(iCol) = cellfun(@(x) replaceValue(x), currentColumn, 'UniformOutput', false);
    elseif isstring(currentColumn) || ischar(currentColumn)
        tab.(iCol) = str2double(currentColumn);
        tab.(iCol)(isnan(tab.(iCol))) = currentColumn(isnan(tab.(iCol)));
    elseif iscategorical(currentColumn)
        continue 
    else
        tab.(iCol)(tab.(iCol) == 555 | tab.(iCol) == 777 | tab.(iCol) == 999) = NaN;
    end
end

% define a function to replace values of "555" with NaN
function out = replaceValue(x)
    if isequal(x, '555') | isequal(x, '777') | isequal(x, '999')
        out = NaN;
    else
        out = x;
    end
end

end 

