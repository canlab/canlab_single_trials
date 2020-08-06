% function [obj, bincnt] = quantileByY(obj, study_id, quantiles)
%
% study_id  - numeric vector of length(obj.Y) elements. Each element
%               contains corresponding subject id (numeric)
% quantiles - Number of quantiles. e.g. 4 means averaging into quartiles.
%
% returns ::
%
% obj       - fmri_data object of size quantiles * length(unique(study_id))
%
% bincnt    - returns number of trials averaged for each quantile. If data
%               is unbalanced, this gives you a rough relative precision
%               estimate for each resultant averaged map.
%
% reads obj.Y, subdivides it into quantiles, and averages elements of obj
% within quantile. Returns new dat object that results from this averaging
%
% Quantiles are computed and averaged within subject.

function [obj, bincnt] = quantileByY(obj, subject_id, parts)
    if length(obj.Y) ~= length(subject_id), error('subject_id length mismatch'); end
    uniq_sid = unique(subject_id);
    
    newdat = cell(length(uniq_sid),parts); % we'll straighten this out later
    bincnt = deal(zeros(length(uniq_sid),parts));
    dropped = zeros(length(uniq_sid),1); 
    for i = 1:length(uniq_sid)
        this_sid = uniq_sid(i);
        this_idx = ismember(subject_id, this_sid);
        this_dat = obj.get_wh_image(this_idx);
        
        [this_bincnt, bin] = bincounts(this_dat.Y, parts);
        dropped(i) = length(this_dat.Y) - sum(this_bincnt);
        
        bincnt(i,:) = this_bincnt;
        for j = 1:parts
            newdat{i,j} = mean(this_dat.get_wh_image(bin == j));
            warning('off','fmriDataSt:mean'); % we only want to see this once
            
            newdat{i,j}.images_per_session = mean(this_dat.images_per_session);
            newdat{i,j}.covariate_names = this_dat.covariate_names;
        end
    end
    newdat = reshape(newdat',numel(newdat),1);
    obj = cat(newdat{:});
    obj.images_per_session = [];
    for i = 1:length(newdat)
        obj.images_per_session = [obj.images_per_session; newdat{i}.images_per_session];
    end
    obj.history{end+1} = ['averaged Y over ' int2str(parts) ' quantiles within subject'];
    bincnt = reshape(bincnt',numel(bincnt),1);
    warning('fmriDataSt:mean','on');
end

% quantiles x into # parts.
% [~,bin] = bincounts(x,4)
% Now length(bin) = length(x), with each element indicating if x belongs to
% quartile 1 (lowest), 2, 3 or 4 (highest).
function [bincnt, bin] = bincounts(x,parts)
    % add a bit of random noise to avoid any ties across quartiles
    x = x + rand(length(x),1)*0.0001;

    bincnt = zeros(parts,1);
    bin = zeros(length(x),1);
    breaks = [min(x), quantile(x,parts-1), max(x)];    
    for i = 1:parts
        if i == parts
            idx = x >= breaks(i) & x <= breaks(i+1);
        else
            idx = x >= breaks(i) & x < breaks(i+1);
        end
        bincnt(i) = sum(idx);
        bin(idx) = i;
    end
end
