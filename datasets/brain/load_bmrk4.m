function [image_obj, networknames, imagenames] = load_bmrk4(varargin)
    % This code loads a dataset object saved in a mat file, and attempts to
    % download it if it cannot be found. 
    dataset_name = 'bmrk4';
    
    % we don't use these, but assign them for backwards compatbility
    networknames = {};
    imagenames = {};
    
    checkmd5 = 0;
    verbose = 0;
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'md5check'
                    checkmd5 = 1;
                case 'verbose'
                    verbose = varargin{i+1};
            end
        end
    end
    
    fmri_data_file = which([dataset_name, '_data.mat']);

    if isempty(fmri_data_file)

        fmri_data_file = download_dataset(dataset_name, varargin{:});
        
    end

    if checkmd5 
        if wrongmd5(fmri_data_file,varargin{:})
            fname = dir(fmri_data_file);
            error(['MD5 check failed. Please manually obtain correct dataset, and ensure it is returned',...
                ' correctly when invoking which(''', fname.name, ''')']);
        else
            if verbose, fprintf('MD5 match confirmed.\n'); end
        end
    end
    
    if verbose, fprintf('Loading: %s\n', fmri_data_file); end

    image_obj = fmri_data_st(importdata(fmri_data_file));
    %image_obj.dat = image_obj.dat.*image_obj.metadata_table.globalSignalAmp; % corrects for the failure of single_trial_analysis.m to incorporate global signal scaling information

    if verbose
        descriptives(image_obj); 
        disp('Pain ratings in image_obj.Y');
        disp('Additional metadata in image_obj.additional_info struct');
    end
end
