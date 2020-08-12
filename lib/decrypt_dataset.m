% Front end to facilitate file encryption and decryption of CANlab single 
% trials datasets. Some datasets may be public, but private repositories
% are stored on public webservers in encrypted form. This function decrypts
% them using a presaved decryption key specific to this repository. This
% function uses the Java Cryptox libraries for encryption/decryption.
%
% Input ::
%   inFile  - character array path to encrypted file
%
%   outFile - character array path to location to save output file
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
function decrypt_dataset(inFile, outFile, varargin)

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
    
    cipher.decrypt(inFile, outFile);
end
