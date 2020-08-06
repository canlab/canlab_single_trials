function [obj, objcodes] = cat(obj,varargin)
    % extends fmri_data/cat by merging other fields like metadata_tables
    % as formated in the single trial datasets

    imgcnt = zeros(1 + length(varargin),1);
    imgcnt(1) = size(obj.dat,2);
    for i = 1:length(varargin)
        imgcnt(i+1) = size(varargin{i}.dat,2);
    end
    
    % this notation is necessary to avoid invoking cat recursively
    obj = {obj, varargin{:}};
    
    % fmri_data/cat doesn't concatenate tables at all, so let's do it and
    % save the result for later
    obj = expand_metadata_table(obj);
    metadata_table = cellfun(@(x1)(x1.metadata_table), obj, 'UniformOutput',false);
    
    % fmri_data/cat concatenates additional_info, but this doesn't work if
    % fields aren't all identical, so let's ensure all structs (like
    % additional_info) have identical fields and expand those that lack it
    % with nans
    obj = expand_struct(obj);
    
    [obj, objcodes] = cat@fmri_data(obj{:});
    
    if any(cellfun(@(x1)(~isempty(x1)),metadata_table))
        obj.metadata_table = tbl_vcat(metadata_table{:});
        if size(obj.dat,2) ~= height(obj.metadata_table)
            error('Metdata length mismatch');
        end
    end
    obj.additional_info = merge_struct(obj.additional_info, imgcnt);
end


% returns dat such that all columns of metadata_table match. If one element 
% has metadata_table that lack columns others have we add that column to it 
% and initialize it to a default empty value. This is necessary for table
% concatenation to work.
function dat = expand_metadata_table(dat)
    t = {};
    has_table = zeros(length(dat),1);
    for i = 1:length(dat)
        if ~isempty(dat{i}.metadata_table)
            t{end+1} = dat{i}.metadata_table;
        end
    end

    t = expand_metadata_table(t{:});
    lacks_table = find(~has_table);
    has_table = find(has_table);
    for i = 1:length(t)
        dat{has_table(i)}.metadata_table = t{i};
    end
    varnames = t{1}.Properties.VariableNames;

    numericCol = cellfun(@(x1)(isnumeric(x1)),(table2cell(dat{1}(1,:))))
    charCol = cellfun(@(x1)(isChar(x1)),(table2cell(dat{1}(1,:))))
    if sum(charCol + numericCol) ~= width(dat{1})
        error('Only numeric and character metadata_table entries supported');
    end
    numericCol = find(numericCol); % can assume everything else is character
    for i = 1:length(lacks_table)
        n = size(dat{lacks_table(i)}.dat,2);
        if numericCol(1)
            t = table(nan(n,1),'VariableName',varnames(1));
        else
            t = table(cellstr(repmat('',n,1)),'VariableNames',varnames(1));
        end
        for 2:length(varnames)
            if numericCol(1)
                t = table(nan(size(n,1),'VariableName',varnames(1));
            else
                t = table(cellstr(repmat('',n,1)),'VariableNames',varnames(1));
            end
    end
end



% generate matching fields across all elements of struct-array dat. If a
% field is missing in one element, create it and populate it with a
% default. Needed for struct concatenation to work.
function dat = expand_struct(dat)
    fnames = {};
    isStruct = [];
    for i = 1:length(dat)
        newfnames = fieldnames(dat{i});
        fnames = [fnames, newfnames];
        for j = 1:length(newfnames)
            if isstruct(dat{i}.(newfnames{j}))
                isStruct = [isStruct; 1];
            else
                isStruct = [isStruct; 0];
            end
        end
    end
    [uniq_fnames, b] = unique(fnames);
    isStruct = isStruct(b);
    
    % recurse on any structs
    for i = 1:length(uniq_fnames)
        if isStruct(i)
            new_struct = cell(length(dat),1);
            for j = 1:length(dat)
                if ~ismember(uniq_fnames{i},properties(dat{j}))
                    dat{j}.(uniq_fnames{i}) = struct('');
                end
                new_struct{j} = dat{j}.(uniq_fnames{i});
            end
            new_struct = expand_struct(new_struct);
            for j = 1:length(dat)
                dat{j}.(uniq_fnames{i}) = new_struct{j};
            end
        end
    end
    
    % pad missing elements with nans
    for i = 1:length(dat)
        these_fnames = fieldnames(dat{i});
        missing = find(~ismember(uniq_fnames, these_fnames));
        for j = 1:length(missing)
            this_fname = uniq_fnames(missing(j));
            dat{i}.(this_fname) = nan;
        end
    end
end
    
% concatenates fields of struct-array elements. If one of the fields is 
% also a struct we recurse.
function ad_info = merge_struct(ad_info, imgcnt)
    n = length(ad_info);
    fnames = fieldnames(ad_info);
    for i = 1:length(fnames)
        dat = cell(n,1);
        % pull each element (e.g. all ratings elements, all temperature 
        % elements, or whatever we're dealing with, and load them into a 
        % cell array.
        for j = 1:n 
            dat{j} = ad_info(j).(fnames{i});
        end
        % figure out an appropriate concatenation procedure
        if isa(dat{1},'struct') % recurse
            ad_info(1).(fnames{i}) = merge_struct(dat{:}, imgcnt);
        elseif isa(dat{1},'char') 
            % study_id and subject_id most likely. This will be vertical,
            % because a horizontal cell array is just a string and we want 
            % to keep entries from different constitutent datasets as
            % separate elements.
            ad_info(1).(fnames{i}) = char(dat);
        else % most likely numeric
            if sum(cellfun(@(x1)(size(x1,1)),dat)) == sum(imgcnt)
                catMethod = 'vertical';
            elseif sum(cellfun(@(x1)(size(x1,2)),dat)) == sum(imgcnt)
                catMethod = 'horizontal';
            else
                [nn,mm] = size(dat{1});
                if nn >= mm % vertical concat
                    if nn == mm
                        warning(['''' fnames{i}, ''' is square, cannot autodetermine concatenation orientation. Defaulting to vertical concatenation. Check results.']);
                    end
                    catMethod = 'vertical';
                else
                    if isnumeric(dat{1})
                        warning(['Only tested with tall numeric data but ' fnames{i} ' is not tall.']);
                    end
                    catMethod = 'horizontal';
                end
            end
               
            switch catMethod
                case 'vertical'
                    ad_info(1).(fnames{i}) = cat(1,dat{:});
                case 'horizontal'
                    ad_info(1).(fnames{i}) = cat(2,dat{:});
                otherwise
                    warning('unhandeled datatype in struct. Erasing entry');
                    ad_info(1).(fnames{i}) = [];
            end
        end
    end
    if ~isempty(ad_info)
        ad_info = ad_info(1);
    end
end
