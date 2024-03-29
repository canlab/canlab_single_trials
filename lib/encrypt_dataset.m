% Front end to facilitate file encryption and decryption of CANlab single 
% trials datasets. Some datasets may be public, but private repositories
% are can be stored on public webservers in encrypted form. This function 
% facilitates creation of such encrypted files for use with this
% repository. It will use a presaved key so that decryption integrates 
% seamlessly with the repository. After uploading an encrypted file to a
% public webserver, add an entry in download_dataset to download that file
% and use decrypt_dataset() to decode it.
%
% This function uses the Java Cryptox libraries for encryption/decryption.
%
% Input ::
%   inFile  - character array path to unencrypted file
%
%   outFile - character array path to location to save encrypted file
%
% Optional Input ::
%
%   'key'   - followed by a string specifying the password to use for decryption
%
% Dependencies ::
%
%  AES class
%   - javaMethodWrapper class
%
function encrypt_dataset(inFile, outFile, varargin)

    secret = [];
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'key'
                    secret = varargin{i+1};
            end
        end
    end

    if isempty(secret)
        keyfile = fopen('single_trials_aes_key.txt');
        secret = fscanf(keyfile,'%s');
        fclose(keyfile);
    end
    
    cipher = AES(secret,'SHA-256');
    
    cipher.encrypt(inFile, outFile);
end
