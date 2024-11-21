% overloads fmri_data/mean and adds on functionality to average
% fmri_data_st.Y, *.covariates, *.metadata_table
function [newobj, varargout] = mean(obj, varargin)
    % we need this to avoid image_vector/mean introducing nans
    obj = obj.remove_empty;
    % I don't know how to use varargout, and someone who knows more should
    % replace this switch with a more elegant solution
    switch nargout
        case 0 
            newobj = mean@fmri_data(obj, varargin{:});
        case 1
            newobj = mean@fmri_data(obj, varargin{:});
        case 2
            [newobj, varargout{1}] = mean@fmri_data(obj, varargin{:});
        case 3
            [newobj, varargout{1}, varargout{2}] = mean@fmri_data(obj, varargin{:});
        otherwise
            error('Does not support more than 3 output arguments');
    end
    
% All the below is now integrated into fmri_data.mean.  tor wager, 12/5/2023    
%     copyfield = {'source_notes','images_per_session','Y_names','Y_descrip',...
%         'covariate_names','additional_info','history'};
%     for i = 1:length(copyfield)
%         newobj.(copyfield{i}) = obj.(copyfield{i});
%     end
%     
%     newobj.Y = nanmean(obj.Y);
%     newobj.covariates = nanmean(obj.covariates,2);
%     
%     t = obj.metadata_table;
%     if ~isempty(t)
%         tnames = t.Properties.VariableNames;
%         vars = cell(1,length(tnames));
%         for i = 1:length(tnames)
%             if isnumeric(t.(tnames{i}))
%                 vars{i} = nanmean(t.(tnames{i}));
%             elseif iscell(t.(tnames{i}))
%                 vars{i} = t.(tnames{i})(1);
%             elseif ischar(t.(tnames{i}))
%                 vars{i} = t.(tnames{i})(1,:);
%             else
%                 warning(['No policy for averaging datatype ''' class(t.(tnames{i})) ''' in metadata_table. Dropping column ''' tnames{i} '''']);
%                 vars{i} = [];
%             end
%         end
%         newobj.metadata_table = cell2table(vars,'VariableNames',tnames);
%         %warning('fmriDataSt:mean','Any metadata_table char arrays converted to cellstr');
%     end
%     
%     % mean@fmri_data recasts this, so let's fix that here
%     newobj = fmri_data_st(newobj);
end % function
