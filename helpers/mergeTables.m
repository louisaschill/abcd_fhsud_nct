function mergedTable = mergeTables(tabs)
    % read in the first table to use as a starting point for the merge
    mergedTable = tabs{1};
    
    % loop through the remaining tables and merge them one by one
    for i = 2:length(tabs)
        % read in the next table to merge
        nextTable = tabs{i};
        
        % merge the next table with the current merged table
        mergedTable = outerjoin(mergedTable, nextTable, 'Keys','subjectkey','Type', 'left','MergeKeys', true);
    end
end