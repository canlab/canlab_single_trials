% returns dat such that all columns of metadata_table match. If one element 
% has metadata_table that lack columns others have we add that column to it 
% and initialize it to a default empty value. This is necessary for table
% concatenation to work.
function t = expand_metadata_table(varargin)
    % identify all fields in all structs, determine if lengths match image 
    % count and if so determine the orientation of the matching dimension.
    % also identify fields that are structs for recursing later
    [fnames, fclass] = deal([]);
    for i = 1:length(varargin)
        these_fnames = varargin{i}.Properties.VariableNames;
        fnames = [fnames, these_fnames];
        for j = 1:length(these_fnames)
            fclass = [fclass, {class(varargin{i}.(these_fnames{j}))}];
        end
    end
    [uniq_fnames, b] = unique(fnames);
    uniq_fclass = fclass(b);
    
    if isempty(uniq_fnames) % this is the empty field scenario
        return
    end
    
    % sanity check to ensure all instances of uniq_fnames have the same
    % class
    for i = 1:length(uniq_fnames)
        instances = find(strcmp(uniq_fnames{i},fnames));
        
        if ~all(strcmp(fclass(instances),fclass(instances(1))))
            error([uniq_fnames{i}, ' does not have the same datatype across tables']);
        end
    end
    
    % for each these_structs see if field is missing
    % if missing replace with nan or class equivalent
    t = cell(length(varargin),1);
    for i = 1:length(t)
        t{i} = varargin{i};
        hasCol = t{i}.Properties.VariableNames;
        
        missCol = find(~ismember(uniq_fnames, hasCol));
        for j = 1:length(missCol)
            this_col = uniq_fnames{missCol(j)};
            this_class = uniq_fclass{missCol(j)};
            if strcmp(this_class,'cell')
                val = cell(height(t{i}),1);
            else
                val = nan(height(t{i}),1);
                val = eval([this_class, '(val);']); % cast to approprite datatype
            end
            t2 = array2table(val,'VariableNames',{this_col});
            t{i} = [t{i},t2];
        end
    end
end