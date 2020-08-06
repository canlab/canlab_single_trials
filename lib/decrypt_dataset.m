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
% Dependencies ::
%
%  AES class
%   - javaMethodWrapper class
%
function decrypt_dataset(inFile, outFile)
    keyfile = fopen('single_trials_aes_key.txt');
    secret = fscanf(keyfile,'%s');
    fclose(keyfile);
    
    cipher = AES(secret,'SHA-256');
    
    cipher.decrypt(inFile, outFile);
end