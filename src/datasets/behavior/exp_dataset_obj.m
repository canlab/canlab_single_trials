clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2

fprintf('Preparing EXP dataset obj...\n');
dataset_obj = canlab_dataset;
dataset_obj.Description.Experiment_Name = 'exp';
dataset_obj.Description.references = {'Atlas, et al. (2010) Journal of Neuroscience'};

%% update these paths as needed
dataRoot1     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/Expectancy/';
% single trial image metadata
imgmeta        = importdata([dataRoot1 '/exp_meta.mat']); 
% behavioral meta (might be identical to the above)
meta = importdata('/work/ics/data/projects/wagerlab/labdata/single_trials_experimental/Data/Expectancy/exp_meta.mat');
dataRoot = ['/work/ics/data/projects/wagerlab/labdata/current/Expectancy'];
cue_data = importdata([dataRoot, '/Behavioral/Sas_proc_mixed/behav_data_byrun_subjsusing.mat']);


%% ideally these will match across studies
dataset_obj.Subj_Level.names = {'age','male','race','right_handed'};
dataset_obj.Subj_Level.type = {'int','boolean','level','boolean'};
dataset_obj.Event_Level.names = {'rating', 'T', 'soundintensity', 'cue', ...
    'social', 'placebo', 'value', 'runN', 'siteN', 'overallN', 'vif', ...
    'ctrl', 'reveal','conditioningN','reg','handholding','drug','open'};
dataset_obj.Event_Level.type = ...
    repmat('numeric',length(dataset_obj.Description.Event_Level));
dataset_obj.Event_Level.units = ...
    {'pixels','C','levels','levels','levels','boolean','levels','count','count',...
    'count','unitless','boolean','boolean','count','levels','boolean','boolean','boolean'};

%% some study specific stuff
datset_obj.Event_Level.descrip{ismember(dataset_obj.Event_Level.names,'rating')} = ...
    'VAS is 0-600 here, which was determined by screen resolution (units = pixels)';

%% get within subject data
nSubj = length(imgmeta.dat_obj);
hasData = zeros(nSubj,1);
subjVols = cell(nSubj,1);
volsPerSubj = zeros(nSubj,1);


for i=1:nSubj
    %% load subject's fmri_data obj
    subjNum = i;
    subj = imgmeta.dat_obj{subjNum};
    dataset_obj.Subj_Level.id{end + 1} = strrep(subj,'.mat','');
    %fprintf('Processing subject %s ...\n', subj);

    file = [dataRoot1,subj];
    fileExt = strsplit(file,'.');
    if strcmp(fileExt{end},'mat')
        img = importdata(file);
    else
        if strcmp(fileExt{end},'nii')
            img = fmri_data(file);
        end
    end

    %for levoderm data
    if isfield(img,'dat')
        img = img.dat;
    end

    %{
    %% record number of vols pertaining to this subj
    % comment this out for original nsf data
    if isempty(img.images_per_session)
         img.images_per_session = ones(1,size(img.dat,2));
    end
    %}
    
    subjY = imgmeta.ratings{i};
    if length(subjY) ~= size(img.dat,2)
        error('metadata mismatch');
    end
    subjYIdx = ~isnan(subjY);
    cueIdx = 0;
    for k = 1:length(cue_data.X1data.ratings)
        candidate = cue_data.X1data.ratings{k};
        candidateIdx = ~isnan(candidate);
        if sum(candidateIdx) == sum(subjYIdx)
            r = corr(candidate(candidateIdx),subjY(subjYIdx));
            if r + 0.001 > 1
                cueIdx = k;
            end
        else
            n = min([length(candidateIdx),length(subjYIdx)]);
            if sum(candidateIdx(1:n)) == sum(subjYIdx(1:n))
                r = corr(candidate(candidateIdx(1:n)),subjY(subjYIdx(1:n)));
                if r + 0.001 > 1
                    cueIdx = k;
                end
            end
        end
    end
    these_cov = [];
    if cueIdx > 0
        these_cov = zeros(length(cue_data.X1data.ratings{cueIdx}), length(dataset_obj.Event_Level.names));
        candidateIdx = find(isnan(cue_data.X1data.ratings{cueIdx}) ~= 1);

        X1 = cue_data.X1{cueIdx};
        T = cue_data.X1data.temps{cueIdx};
        % rating, temp, cue, placebo, value, runN
        these_cov(:,ismember(dataset_obj.Event_Level.names,{'rating'})) = cue_data.X1data.ratings{cueIdx};
        these_cov(:,ismember(dataset_obj.Event_Level.names,{'T'})) = T;
        these_cov(:,ismember(dataset_obj.Event_Level.names,{'cue'})) = X1(:,1)/2;
        these_cov(:,ismember(dataset_obj.Event_Level.names,{'runN'})) = X1(:,4);

        seq = [];
        runs = unique(these_cov(:,ismember(dataset_obj.Event_Level.names,{'runN'})));
        for k = 1:length(runs)
            this_run = runs(k);
            run_idx = find(these_cov(:,ismember(dataset_obj.Event_Level.names,{'runN'})) == this_run);
            seq = [seq; (1:length(run_idx))'];
        end        
        these_cov(:,ismember(dataset_obj.Event_Level.names,{'siteN'})) = seq(:);
        these_cov(:,ismember(dataset_obj.Event_Level.names,{'overallN'})) = (1:size(these_cov,1))';
        %these_cov = these_cov(candidateIdx,:);

        meta_idx = ~isnan(meta.ratings{i});
        idx = ~isnan(these_cov(:,1));
        
        if sum(idx) ~= sum(meta_idx)
            n = min([length(idx),length(meta_idx)]);
            if sum(idx(1:n)) == sum(meta_idx(1:n))
                r = corr(these_cov(idx(1:n),1),meta.ratings{i}(meta_idx(1:n)));
                if r + 0.001 > 1
                    idx = idx(1:n);
                    meta_idx = meta_idx(1:n);
                    these_cov = these_cov(1:n,:);
                end
            end
        end
        
        vif = zeros(size(these_cov,1),1);
        if corr(meta.ratings{i}(meta_idx),these_cov(idx,1)) + 0.0001 >= 1
            vif(meta.high_vif_trials_idx{i}(idx)) = 2.5;
            these_cov(:,ismember(dataset_obj.Event_Level.names,{'vif'})) = vif;
        else
            warning(['Could not recover VIF for subject ' int2str(i)]);
        end
    else
        warning(['Could not find ratings for subject ' num2str(i)]);
        keyboard
    end
    dataset_obj.Event_Level.data{end+1} = these_cov;
end

% drop empty fields
hasData = any(cell2mat(dataset_obj.Event_Level.data'));
for i = 1:length(dataset_obj.Event_Level.data)
    dataset_obj.Event_Level.data{i} = dataset_obj.Event_Level.data{i}(:,hasData);
end
dataset_obj.Event_Level.names = dataset_obj.Event_Level.names(hasData);
dataset_obj.Event_Level.type = dataset_obj.Event_Level.type(hasData);
dataset_obj.Event_Level.units = dataset_obj.Event_Level.units(hasData);

if length(dataset_obj.Event_Level.data) ~= nSubj
    error('Length of dataset mismatch');
else
    save([analysisRoot 'resources/prep_canlab_dataset_objs/' ...
        dataset_obj.Description.Experiment_Name, '_dataset_obj.mat'],'dataset_obj');
end
% perform some kind of sanity check (e.g. positive temperature effects,
% effective experimental manipulations, etc.) just to be safe

covariates = cell2mat(dataset_obj.Event_Level.data');
% add subj id column
sids = [];
for i = 1:length(dataset_obj.Event_Level.data)
    sids = [sids; i*ones(size(dataset_obj.Event_Level.data{i},1),1)];
end
covariates = [covariates, sids];
covariate_names = [dataset_obj.Event_Level.names, 'sid'];

% center within subject
% T - we only care about increments over pain threshold, so units C arent
%   helpful. Better to center
% rating - centered because scale is arbitrary, and this saves us from
%   needing to model intercept
% siteN - centered because we care about effects in the middle trial, not
%   at the zeroth trial
% overallN - ditto
uniq_sids = unique(sids);
cmat = covariates;
for i = 1:length(uniq_sids)
    this_sid = uniq_sids(i);
    pain_idx = find(this_sid == sids);
    mean_v = mean(covariates(pain_idx,:));
    column_idx = find(ismember(covariate_names,{'rating','T','siteN','overallN'}));
    cmat(pain_idx,column_idx) = cmat(pain_idx,column_idx) - repmat(mean_v(column_idx),length(pain_idx),1);
end

% sanity check for order and cue effects
tbl = array2table(cmat,'VariableNames',covariate_names);
m = fitlme(tbl,'rating ~ T + cue + siteN - 1 + (T + cue +siteN - 1 | sid) + (T + cue + siteN | runN)','FitMethod','REML')

%{
m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations             887
    Fixed effects coefficients           3
    Random effects coefficients         74
    Covariance parameters               17

Formula:
    rating ~ T + cue + siteN + (T + cue + siteN | sid) + (1 + T + cue + siteN | runN)

Model fit statistics:
    AIC       BIC       LogLikelihood    Deviance
    9777.6    9873.3    -4868.8          9737.6  

Fixed effects coefficients (95% CIs):
    Name           Estimate    SE        tStat      DF     pValue        Lower      Upper  
    'T'             27.997     2.9009     9.6514    884             0     22.304     33.691
    'cue'           106.17     13.512     7.8574    884    1.1324e-14     79.649     132.69
    'siteN'        -9.3379     1.0648    -8.7699    884             0    -11.428    -7.2481

Random effects covariance parameters (95% CIs):
Group: sid (14 Levels)
    Name1          Name2          Type          Estimate    Lower       Upper  
    'T'            'T'            'std'           8.5042       4.909     14.732
    'cue'          'T'            'corr'        0.017015    -0.62671    0.64694
    'siteN'        'T'            'corr'        -0.63314    -0.99154    0.84475
    'cue'          'cue'          'std'           39.475      23.572     66.108
    'siteN'        'cue'          'corr'        -0.75572    -0.99868    0.93423
    'siteN'        'siteN'        'std'           1.9147     0.45413     8.0725

Group: runN (8 Levels)
    Name1                Name2                Type          Estimate    Lower       Upper   
    '(Intercept)'        '(Intercept)'        'std'           8.1905      4.0791      16.446
    'T'                  '(Intercept)'        'corr'        -0.80196         NaN         NaN
    'cue'                '(Intercept)'        'corr'          0.7077     0.70458     0.71079
    'siteN'              '(Intercept)'        'corr'        -0.92082    -0.92847    -0.91239
    'T'                  'T'                  'std'           4.2451      1.7711      10.175
    'cue'                'T'                  'corr'         -0.9896     -0.9902    -0.98897
    'siteN'              'T'                  'corr'         0.97143     0.96967     0.97308
    'cue'                'cue'                'std'           18.413      6.1496       55.13
    'siteN'              'cue'                'corr'        -0.92719    -0.93192    -0.92215
    'siteN'              'siteN'              'std'           1.7903     0.60966      5.2576

Group: Error
    Name             Estimate    Lower     Upper 
    'Res Std'        56.877      54.207    59.678
%}