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
        if isempty(which([dataset_name, '_data.mat']))
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
    if verbose, fprintf('Downloading %s...\n', outFile); end
    switch dataset_name
        case 'nsf'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22472279',...
                'private_link','aeb02d8d0653bce75a48');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'bmrk3pain'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22472291',...
                'private_link','eac57f22a5c6b57b76f8');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'bmrk3warm'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22472312',...
                'private_link','4d5bde58f074b33e2087');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'bmrk4'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22472405',...
                'private_link','6844caa3143cc3f8065e');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'bmrk5pain'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22483103',...
                'private_link','005d207aa94f57c09326');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'bmrk5snd'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22472615',...
                'private_link','372e18b2928ec9c93291');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'remi'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22473539',...
                'private_link','a7ee52241d08eccf6a14');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'scebl'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22473647',...
                'private_link','8fe187af438304411d35');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'ie2'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22473053',...
                'private_link','3ba42d32113a98ed0252');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'ie'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22472924',...
                'private_link','dcf2fd298d625526c8f6');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'exp'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22473785',...
                'private_link','da2c2f8ad2635e489d5d');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'levoderm'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22473455',...
                'private_link','7eb4ef2a9cd0a53fd5d3');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'stephan'
            % this is wrong, update after upload
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22475531',...
                'private_link','ce68f0116425bfca06a1');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'romantic'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22473578',...
                'private_link','3539b20d01cf0c342bd5');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        case 'ilcp'
            inFile = websave('tmp.mat_encrypted','https://ndownloader.figshare.com/files/22473332',...
                'private_link','36c275944d646163e4b1');
            decrypt_dataset(inFile, outFile);
            delete tmp.mat_encrypted
        otherwise
            error(['Dataset ', dataset_name, ' download not configured']);
    end
    
    path = which(outFile);
end