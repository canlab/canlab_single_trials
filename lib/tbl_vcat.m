% table/vertcat doesn't handle string concatenation well, so we use this
% instead. It assumes tables already all have the same columns.
% this should probably be rewritten as a recursive function instead of a 
% varargin function
%
% Input ::
%
%   vararing    - variable length list of tables. e.g. tbl_vcat(tbl1,
%                   tbl2,..., tbln)
%
% Written by Bogdan Petre, April 2020
function t = tbl_vcat(varargin)
    fnames = {};
    for i = 1:length(varargin)
        fnames = [fnames, varargin{i}.Properties.VariableNames];
    end
    fnames = unique(fnames);
    for i = 1:length(varargin)
        hasfnames = ismember(fnames, varargin{i}.Properties.VariableNames);
        if ~all(hasfnames)
            error('Tables must have identical column names for vertcat()');
        end
    end
    
    colVal = cell(1,length(fnames));
    for i = 1:length(fnames)
        fclass = class(varargin{1}.(fnames{i}));
        colVal{i} = cellfun(@(x1)(x1.(fnames{i})),varargin,'Uniform',false)';
        
        if ischar(colVal{i}{1})
            colVal{i} = char(colVal{i}{:});
        else
            colVal{i} = cat(1,colVal{i}{:});
        end
    end
    t = table(colVal{:},'VariableNames',fnames);
end