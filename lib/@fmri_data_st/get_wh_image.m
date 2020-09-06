% function dat = get_wh_image(dat, idx)
%
% overloading of internal canlab fxn, because canlab function doesn't
% handle all fields of fmridat_obj correctly with it's own get_wh_image
function obj = get_wh_image(obj, idx)
    if iscell(obj.additional_info)
        obj.additional_info = cell2mat(obj.additional_info');
    end
    imgcnt = size(obj.dat,2);
    obj = get_wh_image@fmri_data(obj,idx);
    
    % this is the stuff the canlab function does not do
    obj = slice_struct(obj, imgcnt, idx);
end

% a recursive function that slices all unsliced elements of the image data,
% (a struct), but if one of its elements is itself a struct then we recurse
% and apply the same logic to it.
function newdat = slice_struct(dat, imgcnt, idx, varargin)
    newdat = dat;
    fnames = fieldnames(dat);
    top_fnames = fieldnames(fmri_data)';
    for i = 1:length(fnames)
        this_fname = fnames{i};
        this_field = dat.(this_fname);
        if isstruct(this_field)
            if ~ismember(this_fname,{'mask','volInfo'})
                ffnames = {varargin{:}, this_fname}; % this is used for throwing informative errors
                newdat.(this_fname) = slice_struct(this_field, imgcnt, idx, ffnames{:});
            end
        elseif ~ismember(this_fname,top_fnames) % prevent it from modifying toplvl fields
            ffnames = {varargin{:}, this_fname}; % this is needed for throwing informative errors
            newdat.(this_fname) = get_slice(this_field, imgcnt, idx, ffnames{:});
        end
    end
end

function slice = get_slice(dat, imgcnt, idx, varargin)
    [n,m] = size(dat);
    
    if n == imgcnt
        % this prints out the field and any parent fields of ambiguous slicing
        if n == m
            errMsg = ['Connot determine orientation of '];
            for i = 1:length(varargin)
                errMsg = [errMsg, varargin{i}, ' > '];
            end
            warning([errMsg, '. Assuming rows are trials.\n']);
        end

        slice = dat(idx,:);
    elseif m == imgcnt
        slice = dat(:,idx);
    else
        slice = dat;
    end
end