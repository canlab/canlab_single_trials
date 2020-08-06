classdef AES < handle
    % AES implementation of file encryption/decryption with prespecified key
    %   AES is an industry standard cipher for password protecting data.
    %   This class provides a way of initializing a particular implementation
    %   of this cipher optimized for encryption of large files (although it
    %   should also work for smaller files).
    %
    % A password hash must be provided. If reusing this password in multiple
    %   locations you should rewrite portions of this class to be invoked 
    %   with the password hash itself, not the plaintext password, but for
    %   the purposes of this repo the password was purpose made and is not
    %   reused elsewhere, so for simplicity it's been stored in plain text.
    %   It wouldn't be too hard to modify to take a hash as input instead of
    %   the password itself.
    %
    % Available hashing algorithms are listed here I believe,
    %   https://docs.oracle.com/javase/8/docs/technotes/guides/security/StandardNames.html#MessageDigest
    %
    % Example ::
    %
    %   **Encryption**
    %
    %   coder = AES('thisismypassword', 'SHA-256');
    %   inFile = '/home/bogdan/myLargeFile.bin';
    %   outFile = '/home/bogdan/myLargeFile.bin_encrypted';
    %   AES.encrypt(inFile, outFile);
    %
    %  **Decryption**
    %  
    %  coder = AES('thisismypassword', 'SHA-256');
    %  inFile = '/home/bogdan/myLargeFile.bin_encrypted'
    %  outfile = '/tmp/myLargeFile.bin'
    %  AES.decrypt(inFile, outFile);
    %
    %  % Now /tmp/myLargeFile.bin should be identical to
    %  % /home/bogdan/myLargeFile.bin and you should be able to confirm this
    %  % with an md5sum check on the command line (or in the case of this
    %  % canlab_single_trials repo whic this implementaiton was designed for
    %  % you could use the GetMD5 matlab function).
    %
    % Dependencies ::
    %
    %  - javaMethodWrapper class (available in canlab_single_trials repo)
    %
    % Writen by Bogdan Petre, April, 2020, building off of 
    % https://www.mathworks.com/matlabcentral/fileexchange/73037-matlab-aes-encryption-decryption-example
    
    properties (Access = private)
        secretKey
        cipher
    end
    
    methods
        function obj = AES(secret, algorithm)
            %AES Construct an instance of this class
            %   algorithm options are https://docs.oracle.com/javase/9/docs/specs/security/standard-names.html#messagedigest-algorithms
            import java.security.MessageDigest;
            import java.lang.String;
            import java.util.Arrays;
            import javax.crypto.Cipher;
            
            key = String(secret).getBytes("UTF-8");
            sha = MessageDigest.getInstance(algorithm);
            key = sha.digest(key);
            key = Arrays.copyOf(key, 16);
            obj.secretKey = javaObject('javax.crypto.spec.SecretKeySpec',key, "AES");
            obj.cipher = Cipher.getInstance("AES/CBC/PKCS5Padding");
            %obj.cipher = Cipher.getInstance("AES/CBC/NoPadding");
        end
        
        % Written by Bogdan
        function encrypt(obj, inFilePath, outFilePath)
            import java.util.Base64;
            import javax.crypto.Cipher;
            import java.io.File;
            import java.io.FileInputStream;
            import java.io.FileOutputStream;
            import javax.crypto.spec.IvParameterSpec;
            import java.security.SecureRandom;

            % the use of JavaMethodWrapper here is a hack because we can't
            % normally invoked pointers to primative java datatypes in 
            % matlab, but we need these pointers for the I/O buffer to work
            % in our copy loop below. This solution is courtesy of Benjamin 
            % Davis who wrote the JavaMethodWrapper class
            
            % What this achieves is to allow us to call something like
            % byte[] byteSream = new byte[(int) someLength];
            % inputStream.read(byteStream)
            inputFile = File(inFilePath);
            inputStream  = FileInputStream(inputFile);
            inputMethodObj = JavaMethodWrapper(inputStream, 'read(byte[])');

            % What this achieves is to allow us to call something like
            % outputStream.write(byteStream, 0, decrypted.length())
            outputFile = File(outFilePath);
            outputStream = FileOutputStream(outputFile);
            outputMethodObj = JavaMethodWrapper(outputStream, 'write(byte[], int, int)');

            % get random initialization vector (IV), using similar
            % principles as above
            ivGen = SecureRandom();
            ivGenMethodObj = JavaMethodWrapper(ivGen, 'nextBytes(byte[])');
            iv = zeros(1,16,'int8');
            [~,iv] = ivGenMethodObj.invoke(ivGen,iv);
            
            % lets write the IV to the start of our output file
            outputMethodObj.invoke(outputStream, iv, int32(0), int32(length(iv)));

            % now let's invoke a copy loop. This pulls data 2^26 bytes at a 
            % time, encryptes them, and copies them to the output file. The 
            % size of these files often exceed the Java heap, which is why 
            % we don't just pull all of the file at once.
            byteStream = zeros(1, 2^26, 'int8');                                                % init 64MB buffer
            [len, byteStream] = inputMethodObj.invoke(inputStream, byteStream);                % read
            byteStream = byteStream(1:len);
            while len > -1
                ivPSpec = IvParameterSpec(iv);                                                  % note use of IV
                obj.cipher.init(Cipher.ENCRYPT_MODE, obj.secretKey, ivPSpec);                   % init cipher
                
                encrypted = obj.cipher.doFinal(byteStream);                                     % encrypt
                outputMethodObj.invoke(outputStream, encrypted, int32(0), int32(length(encrypted)));    % write

                % use last 128 bits (16 bytes) of last byteStream as IV for next key
                iv = byteStream(end-15:end);

                byteStream = zeros(1, 2^26, 'int8');                                            % clear 64MB buffer
                [len, byteStream] = inputMethodObj.invoke(inputStream, byteStream);            % read
                byteStream = byteStream(1:len);
            end
        end
        
        % Written by Bogdan
        function decrypt(obj, inFilePath, outFilePath)
            import java.util.Base64;
            import javax.crypto.Cipher;
            import java.io.File;
            import java.io.FileInputStream;
            import java.io.FileOutputStream;
            import javax.crypto.spec.IvParameterSpec;
            import java.security.SecureRandom;

            % the use of JavaMethodWrapper here is a hack because we can't
            % normally invoked pointers to primative java datatypes in 
            % matlab, but we need these pointers for the I/O buffer to work
            % in our copy loop below. This solution is courtesy of Benjamin 
            % Davis who wrote the JavaMethodWrapper class
            %
            % What this achieves is to allow us to call something like
            % byte[] byteSream = new byte[(int) 22369644];
            % inputStream.read(byteStream)
            inputFile = File(inFilePath);
            inputStream  = FileInputStream(inputFile);
            inputMethodObj = JavaMethodWrapper(inputStream, 'read(byte[])');

            % What this achieves is to allow us to call something like
            % outputStream.read(byteStream, 0, decrypted.length())
            outputFile = File(outFilePath);
            outputStream = FileOutputStream(outputFile);
            outputMethodObj = JavaMethodWrapper(outputStream, 'write(byte[], int, int)');

            % let's get the IV
            iv = zeros(1, 16, 'int8'); % our buffer
            [~, iv] = inputMethodObj.invoke(inputStream, iv);
                                    
            % now let's invoke a copy loop. This pulls data 2^26 bytes 
            % at a time, decrypts them, and copies them to an output file. 
            % The size of these files often exceed the Java heap, so we 
            % can't do it all at once. Hence the loop.
            byteStream = zeros(1, 2^26 + 16, 'int8');                                               % init buffer
            [len, byteStream] = inputMethodObj.invoke(inputStream, byteStream);                    % read
            byteStream = byteStream(1:len);
            while len > -1
                ivPSpec = IvParameterSpec(iv);
                obj.cipher.init(Cipher.DECRYPT_MODE, obj.secretKey, ivPSpec);                       % init cipher with current IV
                     
                decrypted = obj.cipher.doFinal(byteStream);                                 % decrypt
                outputMethodObj.invoke(outputStream, decrypted, int32(0), int32(length(decrypted)));% write

                % get new IV from the tail of what we just decrypted
                iv = decrypted(end-15:end); 
                
                byteStream = zeros(1, 2^26 + 16, 'int8');                                           % clear buffer
                [len, byteStream] = inputMethodObj.invoke(inputStream, byteStream);                % read
                byteStream = byteStream(1:len);
            end
        end
    end
end

