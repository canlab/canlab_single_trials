% checks md5 sum of 'file' against records in checksums.md5 (which should
% be at the same level of the directory tree as this function). Returns
% true if md5 mismatch is found. Otherwise returns false.
%
% Dependencies ::
% 
%   - GetMD5

function check = wrongmd5(file, varargin)
    verbose = 0;
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'verbose'
                    verbose = varargin{i+1};
            end
        end
    end
                    
    if verbose, disp('Performing MD5 check. This may take several minutes for large datasets...'); end
    
    % import all md5 checkums
    fid = fopen('checksums.md5');
    md5s = textscan(fid,'%s');
    fclose(fid);
    
    % identify the md5 sum we want for the current file
    fname = dir(file);
    md5 = md5s{1}{find(contains(md5s{1},fname.name))-1};
    
    % perform the comparison
    x1 = GetMD5(file,'file');
    check = ~strcmp(x1,md5);
end