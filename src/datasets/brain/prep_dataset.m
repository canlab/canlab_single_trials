function [data, hasData] = prep_dataset(dataRoot, metaFileName, varargin)
    %stimulus descriptions are in this meta file
    meta = importdata([dataRoot '/' metaFileName]);

    nSubj = length(meta.dat_obj);
    hasData = zeros(nSubj,1);
    subjVols = cell(nSubj,1);
    
    % canlab_dataset object was provided we'll use it to populate
    % additional_info fields
    dataset_obj = [];
    if ~isempty(varargin)
        dataset_obj = varargin{1};
    end
    %can be parfor, but mind the memory usage
    for i=1:nSubj
        %% load subject's fmri_data obj
        subjNum = i;
        subj = meta.dat_obj{subjNum};

        file = [dataRoot,subj];
        fileExt = strsplit(file,'.');
        if strcmp(fileExt{end},'mat')
            img = importdata(file);
        else
            if strcmp(fileExt{end},'nii')
                img = fmri_data(file);
            end
        end
        
        %for levoderm data
        if isfield(img,'dat')
            img = img.dat;
        end
        
        if ~isempty(dataset_obj)
            ad_info = dataset_obj.Event_Level.data{i}; % event level metadata (e.g. cues, stimulus intensities, etc)
            
            % include all ad_info that have entries (some may be empty
            % columns depending on canlab_dataset object implementation).
            has_ad_info = any(cell2mat(dataset_obj.Event_Level.data'));
            try
                t1 = array2table(ad_info(:,has_ad_info),...
                    'VariableNames', dataset_obj.Event_Level.names(has_ad_info));
            catch
                keyboard
            end
            % add subject_id info
            t2 = table(cellstr(repmat(strrep(subj,'.mat',''), size(img.dat,2),1)),'VariableNames',{'subject_id'});
            t = horzcat(t1,t2);
            img.metadata_table = t;
            
            img.Y = dataset_obj.Event_Level.data{i}(:, ismember(dataset_obj.Event_Level.names,'rating'));
        else
            % none of the original datasets use this statement, it's here 
            %  for future datasets which may lack canlab_dataset objects.
            if size(img.dat,2) == length(meta.ratings{i})
                img.Y = meta.ratings{i};
                img.metadata_table = array2table(...
                    [meta.rating{i}(:), meta.temp{i}(:), repmat(strrep(subj,'.mat',''), size(img.dat,2),1)],...
                    'VariableNames',{'rating','T','subject_id'});
            else
                warning('Could not pair subject with metadata')
                img.Y = nan(size(img.dat,2),1);
                img.metadata_table = array2table(repmat(strrep(subj,'.mat',''), size(img.dat,2),1),...
                    'VariableNames',{'subject_id'});
            end
        end
        % cast for better handling of additional_info field when
        % concatenating
        subjVols{i} = fmri_data_st(img);
    end
    data = cat(subjVols{:});
    
    % do a sanity check that metadata and image lengths match
    if ~isempty(dataset_obj)
        if height(data.metadata_table) ~= size(data.dat,2)
            error('metadata_table dimensions do not match expected number of retained volumes');
        end
        
        if isempty(data.additional_info)
            data.additional_info = struct('references',{dataset_obj.Description.references});
        else
            data.additional_info.references = dataset_obj.Description.references;
        end
    end
    
    % fix some fmri_data object metadata bugs
    data.removed_images = zeros(size(data.dat,2),1);
    % compress if possible
    data = data.remove_empty;
    % cast back to common datatype
    data = fmri_data(data);
end