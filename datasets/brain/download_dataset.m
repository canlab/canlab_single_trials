% Downloads preconfigured datasets from prespecified locations, and
% decrypts them if necessary.
%
% Input ::
%
%   dataset_name    - Name of the dataset to download. Options are: nsf,
%                       bmrk3pain, bmrk3warm, bmrk4, bmrk5pain, bmrk5snd,
%                       remi, scebl, ie2, ie, exp, levoderm, stephan,
%                       romantic, ilcp. Pass as character arrays.
%
% Optional ::
%   
%   'forcedl'         - supress user prompt. Useful for downloading multiple
%                       datasets.
%
%   'verbose'       - followed by 0/1 flag indicating whether to print out
%                       informative text. Default = true
%
% Written by Bogdan Petre 4/29/2020

function path = download_dataset(dataset_name, varargin)

    path = [];
    forcedl = 0;
    verbose = 1;
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'verbose'
                    verbose = varargin{i+1};
                case 'forcedl'
                    forcedl = 1;
            end
        end
    end

    % Prompt the user to continue to avoid accidentally downloading a huge
    % amount of data unknowingly
    if ~forcedl
        if isempty(which([dataset_name, '_data.mat'])) and isempty(which([dataset_name, '_data.nii.gz']))
            fprintf('%s dataset not found in Matlab path.\n', dataset_name);        
            x = input('Would you like to attempt to download it to the current working directory? y/n\n','s');
        else
            fprintf('A %s dataset was found at %s.\n', dataset_name, which([dataset_name, '_data.mat']));
            x = input('Would you like to attempt to download a different copy to the current working directory? y/n\n','s');
        end
        
        if ~ismember(x,{'y','Y','yes','Yes','YES'})
            return
        end
    end
    
    
    outFile = [dataset_name, '_data.mat'];
    if verbose, fprintf('Downloading %s...\n', dataset_name); end
    switch dataset_name
        case 'nsf'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24165545');
        case 'bmrk3pain'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24211439');
        case 'bmrk3warm'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24211832');
        case 'bmrk4'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24212039');
        case 'bmrk5pain'
            inFile = websave(outFile, 'https://figshare.com/ndownloader/files/34638242');
        case 'bmrk5snd'
            inFile = websave(outFile, 'https://figshare.com/ndownloader/files/34638929');
        case 'scebl'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24213212');
        case 'ie2'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24214952');
        case 'ie'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24215003');
        case 'exp'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24212432');
        case 'stephan'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24211892');
        case 'romantic'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24214487');
        case 'ilcp'
            inFile = websave(outFile, 'https://ndownloader.figshare.com/files/24214928');
        case 'baliki_nac_neuron'
            inFile = websave([dataset_name, '_data.nii.gz'], 'https://figshare.com/ndownloader/files/39477889');
        otherwise
            if exist('download_private_dataset')
                % this function is provided by canlab_single_trials_private for in house use
                download_private_dataset(dataset_name,'forcedl','verbose',0);
            else
                error(['Dataset ', dataset_name, ' download not configured']);
            end
    end
    
    path = which(outFile);
end
