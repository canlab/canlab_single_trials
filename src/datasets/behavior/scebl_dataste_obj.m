clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2

fprintf('Preparing SCEBL dataset obj...\n');
dataset_obj = canlab_dataset;
dataset_obj.Description.Experiment_Name = 'scebl';
dataset_obj.Description.references = {'Koban, et al. (2019) Nature Communications'};

%% update these paths as needed
dataRoot     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/SCEBL_single_trial_Leonie/';
% single trial image metadata
imgmeta = importdata([dataRoot,'SCEBLdata_forTor_N26.mat']);

% behavioral meta (might be identical to the above)
evmeta1 = importdata('/work/ics/data/projects/wagerlab/labdata/single_trials_experimental/Data/SCEBL_single_trial_Leonie/SCEBLdata_forTor_N26.mat');
evmeta2 = importdata('/projects/bope9760/pain_roi_prediction_rnd2/resources/SCEBLmri_data.mat');
vif_dat = importdata('/work/ics/data/projects/wagerlab/labdata/single_trials_experimental/Data/SCEBL_single_trial_Leonie/painvifs.mat');

%% ideally these will match across studies
dataset_obj.Subj_Level.names = {'age','male','race','right_handed'};
dataset_obj.Subj_Level.type = {'int','boolean','level','boolean'};
dataset_obj.Event_Level.names = {'rating', 'T', 'soundintensity', 'cue', ...
    'social', 'placebo', 'value', 'runN', 'siteN', 'overallN', 'vif', ...
    'ctrl', 'reveal','conditioningN','reg','handholding','drug','open'};
dataset_obj.Event_Level.type = ...
    repmat('numeric',length(dataset_obj.Description.Event_Level));
dataset_obj.Event_Level.units = ...
    {'vas','C','levels','levels','levels','boolean','levels','count','count',...
    'count','unitless','boolean','boolean','count','levels','boolean','boolean','boolean'};

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

    file = [dataRoot,subj];
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

    subj_idx = find(~isnan(imgmeta.ratings{i}));
    
    name = evmeta1.subjects{i};
    m2_idx = find(ismember(evmeta2.subjname, name));
    trial_idx = find(~isnan(evmeta2.pain{m2_idx}));

    these_cov = [];
    if corr(evmeta2.pain{m2_idx}(trial_idx), imgmeta.ratings{i}(subj_idx)) + 0.001 > 1 && length(imgmeta.ratings{i}) == size(img.dat,2)
        these_cov = zeros(length(evmeta2.pain{m2_idx}), length(dataset_obj.Event_Level.names));
        these_cov(:,ismember(dataset_obj.Event_Level.names,'rating')) = evmeta2.pain{m2_idx};
        these_cov(:,ismember(dataset_obj.Event_Level.names,'T')) = evmeta2.temp{m2_idx};
        these_cov(:,ismember(dataset_obj.Event_Level.names,'cue')) = evmeta2.cues{m2_idx}/2;
        these_cov(:,ismember(dataset_obj.Event_Level.names,'social')) = evmeta2.social{m2_idx}/2;
        these_cov(:,ismember(dataset_obj.Event_Level.names,'runN')) = evmeta2.block{m2_idx};
        
        blocks = unique(evmeta2.block{m2_idx});
        siteN = [];
        for k = 1:length(blocks)
            this_block = blocks(k);
            nn = sum(evmeta2.block{m2_idx} == this_block);
            siteN = [siteN; (1:nn)'];
        end
        these_cov(:,ismember(dataset_obj.Event_Level.names,'siteN')) = siteN;
        these_cov(:,ismember(dataset_obj.Event_Level.names,'overallN')) = mean(1:length(evmeta2.temp{m2_idx}));
        these_cov(:,ismember(dataset_obj.Event_Level.names,'vif')) = vif_dat{i};
    else
        error(['Could not process ' str2num(i)]);
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
sids = unique(covariates(:,ismember(covariate_names,'sid')));
cmat = covariates;
for i = 1:length(sids)
    this_sid = sids(i);
    pain_idx = find(this_sid == covariates(:,ismember(covariate_names,'sid')));
    mean_v = nanmean(covariates(pain_idx,:));
    column_idx = find(ismember(covariate_names,{'rating','T','siteN','overallN'}));
    cmat(pain_idx,column_idx) = cmat(pain_idx,column_idx) - repmat(mean_v(column_idx),length(pain_idx),1);
end

tbl = array2table(cmat,'VariableNames',covariate_names);
m = fitlme(tbl,'rating ~ T + siteN + cue + social - 1 + (T + siteN + cue + social - 1| sid) + (T + siteN | runN)','FitMethod','REML')
%{
m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations            2442
    Fixed effects coefficients           4
    Random effects coefficients        122
    Covariance parameters               17

Formula:
    rating ~ T + cue + social + siteN + (T + cue + social + siteN | sid) + (1 + T + siteN | runN)

Model fit statistics:
    AIC      BIC      LogLikelihood    Deviance
    19551    19673    -9754.7          19509   

Fixed effects coefficients (95% CIs):
    Name            Estimate    SE         tStat      DF      pValue        Lower       Upper    
    'T'               6.6045     1.0205     6.4715    2438    1.1687e-10      4.6033       8.6057
    'cue'             1.3922    0.76491     1.8201    2438      0.068861     -0.1077       2.8922
    'social'          12.407      2.257      5.497    2438    4.2646e-08      7.9809       16.833
    'siteN'         -0.26923    0.12087    -2.2275    2438      0.026005    -0.50625    -0.032219

Random effects covariance parameters (95% CIs):
Group: sid (26 Levels)
    Name1           Name2           Type          Estimate     Lower       Upper   
    'T'             'T'             'std'            4.4556      3.1009      6.4022
    'cue'           'T'             'corr'          0.97134      0.9677     0.97457
    'social'        'T'             'corr'         -0.18525    -0.37609    0.020667
    'siteN'         'T'             'corr'         0.026446    -0.13297     0.18453
    'cue'           'cue'           'std'            1.1928     0.33725      4.2185
    'social'        'cue'           'corr'        -0.052885    -0.25567     0.15436
    'siteN'         'cue'           'corr'         -0.14369    -0.29553    0.015223
    'social'        'social'        'std'            11.203      8.3822      14.974
    'siteN'         'social'        'corr'          0.19243         NaN         NaN
    'siteN'         'siteN'         'std'           0.39205     0.25622     0.59991

Group: runN (6 Levels)
    Name1                Name2                Type          Estimate    Lower       Upper    
    '(Intercept)'        '(Intercept)'        'std'           2.5764      1.4433        4.599
    'T'                  '(Intercept)'        'corr'          0.9547     0.94795      0.96059
    'siteN'              '(Intercept)'        'corr'        -0.45091     -0.5281     -0.36633
    'T'                  'T'                  'std'          0.73279     0.22762       2.3591
    'siteN'              'T'                  'corr'        -0.16488    -0.26186    -0.064599
    'siteN'              'siteN'              'std'          0.20242    0.083517      0.49058

Group: Error
    Name             Estimate    Lower     Upper 
    'Res Std'        12.68       12.323    13.048
%}