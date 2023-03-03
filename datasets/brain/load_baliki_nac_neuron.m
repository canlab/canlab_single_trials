function [image_obj, networknames, imagenames] = load_baliki_nac_neuron(varargin)
    % This code loads a dataset object saved in a mat file, and attempts to
    % download it if it cannot be found. 
    dataset_name = 'baliki_nac_neuron';
    
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
    
    fmri_data_file = which([dataset_name, '_data.nii.gz']);
    metadata_file = which([dataset_name, '_data.csv']);

    if isempty(fmri_data_file)
        fmri_data_file = download_private_dataset(dataset_name, varargin{:});
    end

    if isempty(metadata_file)
        fmri_data_file = download_private_dataset([dataset_name '_metadata'], varargin{:});
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

    image_obj = fmri_data_st(fmri_data_file);
    if exist(fmri_data_file,'file') == 2 && exist(strrep(fmri_data_file,'.nii.gz','.nii'),'file') == 2
        system(sprintf('rm -f %s',strrep(fmri_data_file,'.nii.gz','.nii')));
    end

    image_obj.metadata_table = readtable(metadata_file);
    image_obj.Y = image_obj.metadata_table.ratingPeak;
    image_obj.metadata_table.patient = zeros(height(image_obj.metadata_table),1);
    image_obj.metadata_table.patient(contains(cellstr(image_obj.metadata_table.subject_id), 'cbp')) = 1;
    image_obj.metadata_table.rating = image_obj.metadata_table.ratingPeak;
    
    
    image_obj.additional_info = struct('references', char({'Baliki M, Geha P, Apkarian AV. (2009) "Parsing pain perception between nociceptive representation and magnitude estimation." Journal of Neurophysiology 101(2) 875-87.', ...
					'Baliki M, Geha P, Fields H, Apkarian AV. (2010) "Predicting Value of Pain and Analgesia: Nucleus Accumbens Response to Noxious Stimuli Changes in the Presence of Chronic Pain." Neuron 66(1) 149-160.'}), ...
                    'protocol', char({'Contrasts do not reflect temporal information. Contrast vector AUC is provided in metadata_table. Multiply by contrast maps to convert to BOLD AUC.',...
                'Slice timing information was missing, which was borrowed from the longitudinal subacute backpain study from openpain.org. This was a shot in the dark, ', ...
                'but both protocols were run on the same equipment by the same people and around the same time, so it''s not an unreasonable guess.',...
                'The task inolved online pain ratings, but we do not have the matched visual task available on openpain.org. There are therefore likely visual and motor ',...
                'confounds in the data. ',...
                'Realtime pain ratings are not demarkated by evoking stimulus. Period beginning at stimulus onset and ending 2TR after stimulus offset was used to compute ',...
                'AUC and peak rating.',...
                'Peak rating shows a better relationship with stimulus intensity than AUC, and is the default ''rating'' variable and Y property.'}));
    image_obj.history = {'This data was prepared by B Petre from data retrieved from openpain.org.', ...
				'fmriprep 20.2.3 was used for preprocessing.',...
                'Single trial models were fit with a 1/180 HPF, no timeseries error model, and the statistical design controlled for 6 motion vectors, and CSF. HRF was modeled using SPM defaults.',...
                'Alternative models using AR(1) timeseries errors, 24 motion vectors, or both were considered. The most probable model was selected using random effects bayesian model selection (MACS toolbox).'};

    if verbose
        descriptives(image_obj); 
        disp('Intensity ratings in image_obj.Y');
        disp('Additional metadata in image_obj.metadata_table');
    end
end