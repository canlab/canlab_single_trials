function [t, cov_names] = get_std_tbl_with_fixed_fx(t)
    cov_names = t.Properties.VariableNames;

    % save useful covariate indices and datatypes
    sid_idx = ismember(cov_names,'subject_id');
    covs = table2array(t(:,~sid_idx));
    cov_names = cov_names(~sid_idx);
    sid = t.subject_id; 
    uniq_sid = unique(sid);

    % create subject fixed effects design columns.
    n_sub = length(unique(sid,'rows'));

    Xsub = [];
    for i = 1:n_sub
        if ischar(sid)
            sid = cellstr(sid);
        end
        if iscellstr(sid)
            Xsub = blkdiag(Xsub,ones(sum(strcmp(sid,uniq_sid(i))),1));
        elseif isnumeric(sid)
            Xsub = blkdiag(Xsub,ones(sum(sid == uniq_sid(i)),1));
        else
            error('Only char, cellstr, and numeric sids supported');
        end 
    end 
    sub_labels = cell(1,n_sub);
    for i = 1:n_sub
        sub_labels{i} = ['sub' int2str(i)];
    end

    % standardize everything whle ignoring nans
    covs = [covs, Xsub];
    scovs = nanstd(covs);
    mcovs = nanmean(covs);
    zcovs = (covs - mcovs)./scovs;
    % note that mean is zero for all columns now, so this won't affect
    % regression estimation
    zcovs(isnan(covs)) = 0;
    % rescale to account for increase in size, with addition of nans
    zcovs = zscore(zcovs);
    
    cov_names = [cov_names, sub_labels];

    % get IV std beta coefficients
    t = array2table(zcovs,'VariableNames',cov_names);
end