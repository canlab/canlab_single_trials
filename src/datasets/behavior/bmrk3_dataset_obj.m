clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2

fprintf('Preparing BMRK3 dataset obj...\n');
dataset_obj = canlab_dataset;
dataset_obj.Description.Experiment_Name = 'bmrk3';
dataset_obj.Description.references = {'Wager,  et al. (2013) New England Journal of Medicine', 'Woo et al. (2015) PLoS Biology'};

%% update these paths as needed
dataRoot     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/Tor_bmrk3_datashare/';
% single trial image metadata
imgmeta        = importdata([dataRoot '/bmrk3_single_trial_model.mat']); 
% behavioral meta (might be identical to the above)
evmeta1 = importdata('/work/ics/data/projects/wagerlab/labdata/single_trials_experimental/Data/Tor_bmrk3_datashare/bmrk3_single_trial_model.mat');    
evmeta2 = importdata('/work/ics/data/projects/wagerlab/labdata/current/BMRK3/EXPT.mat');

%% ideally these will match across studies
dataset_obj.Subj_Level.names = {'age','male','race','right_handed'};
dataset_obj.Subj_Level.type = {'int','boolean','level','boolean'};
dataset_obj.Event_Level.names = {'rating', 'T', 'soundintensity', 'cue', ...
    'social', 'placebo', 'value', 'runN', 'siteN', 'overallN', 'high_vif', ...
    'ctrl', 'reveal','conditioningN','reg','handholding','drug','open'};
dataset_obj.Event_Level.type = ...
    repmat('numeric',length(dataset_obj.Description.Event_Level));
dataset_obj.Event_Level.units = ...
    {'vas','C','levels','levels','levels','boolean','levels','count','count',...
    'count','boolean','boolean','boolean','count','levels','boolean','boolean','boolean'};

%% some study specific stuff
datset_obj.Event_Level.descrip{ismember(dataset_obj.Event_Level.names,'high_vif')} = ...
    'vif > 2.5';
datset_obj.Event_Level.descrip{ismember(dataset_obj.Event_Level.names,'rating')} = ...
    '0-100 (thermal); 101-200 (pain)';

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

    % check for mismatches in ev and image data. 
    this_subj = evmeta1.subjects{i};
    these_cov = [];
    meta2_idx = ismember(evmeta2.subjects,this_subj);
    if ~any(meta2_idx) % one of the subjects has a trailing space in meta2 which screws up matching
        meta2_idx = startsWith(evmeta2.subjects,this_subj);
    end
    these_cov = evmeta2.behavior{meta2_idx};
    rating = these_cov(:,ismember(evmeta2.behavior_names,'Rating_painplus100'));
    if any(arrayfun(@(x1,x2)(x1 ~= x2),evmeta1.ratings{i},rating)) % if any ratings don't match
        if length(rating) ~= size(img.dat,2)
            error(['Ratings mismatch for subject ' int2str(i) ' (id = ' this_subj ')']);
        end
    end
    T = these_cov(:,ismember(evmeta2.behavior_names,'Temp'));
    runN = these_cov(:,ismember(evmeta2.behavior_names,'Trial'));
    siteN = these_cov(:,ismember(evmeta2.behavior_names,'Run'));
    overallN = (1:length(rating))';
    sid = i*ones(length(rating),1);
    vif = zeros(length(rating),1);
    vif(evmeta1.high_vif_trials_idx{i}) = 1;
    reg = zeros(length(rating),1);
    reg(these_cov(:,ismember(evmeta2.behavior_names,'ImagineDown')) == 1) = -0.5;
    reg(these_cov(:,ismember(evmeta2.behavior_names,'ImagineUp')) == 1) = 0.5;

    these_cov = zeros(length(rating),length(dataset_obj.Event_Level.names));
    these_cov(:,ismember(dataset_obj.Event_Level.names,'rating')) = rating;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'T')) = T;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'runN')) = runN;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'siteN')) = siteN;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'overallN')) = overallN;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'high_vif')) = vif;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'reg')) = reg;
    
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
covariate_names{contains(covariate_names,'vif')} = 'vif';

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
for i = 1:length(sids) % center 'column_idx' variables (see below) within subject
    this_sid = sids(i);
    pain_idx = find(this_sid == covariates(:,ismember(covariate_names,'sid')));
    mean_v = mean(covariates(pain_idx,:));
    column_idx = find(ismember(covariate_names,{'rating','T','siteN','overallN'}));
    cmat(pain_idx,column_idx) = cmat(pain_idx,column_idx) - repmat(mean_v(column_idx),length(pain_idx),1);
end


%try
%    dat{1}.covariates = covariates;
%    dat{1}.covariate_names = covariate_names;
%end

tbl = array2table(double(cmat),'VariableNames',covariate_names);
m = fitlme(tbl,'rating ~ T + reg + siteN - 1 + (T + reg + siteN - 1| sid) + (T + reg + siteN | runN)','FitMethod','REML')

%{
m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations            3201
    Fixed effects coefficients           3
    Random effects coefficients        143
    Covariance parameters               17

Formula:
    rating ~ T + siteN + reg + (T + siteN + reg | sid) + (1 + T + siteN + reg | runN)

Model fit statistics:
    AIC      BIC      LogLikelihood    Deviance
    31418    31539    -15689           31378   

Fixed effects coefficients (95% CIs):
    Name           Estimate    SE        tStat      DF      pValue        Lower      Upper   
    'T'             22.239     1.2621     17.621    3198             0     19.764      24.713
    'siteN'        -1.3648     0.5651    -2.4151    3198      0.015788    -2.4728    -0.25676
    'reg'           30.402     5.7546     5.2831    3198    1.3552e-07     19.119      41.685

Random effects covariance parameters (95% CIs):
Group: sid (33 Levels)
    Name1          Name2          Type          Estimate     Lower       Upper  
    'T'            'T'            'std'            6.8112      5.1986     8.9242
    'siteN'        'T'            'corr'         -0.26157    -0.58557    0.13448
    'reg'          'T'            'corr'             0.14    -0.26418    0.50235
    'siteN'        'siteN'        'std'            2.8845      2.1482     3.8732
    'reg'          'siteN'        'corr'        -0.099383     -0.4817    0.31471
    'reg'          'reg'          'std'            28.593      20.944     39.036

Group: runN (11 Levels)
    Name1                Name2                Type          Estimate    Lower       Upper   
    '(Intercept)'        '(Intercept)'        'std'           10.164      6.7466      15.312
    'T'                  '(Intercept)'        'corr'        -0.96357    -0.96674     -0.9601
    'siteN'              '(Intercept)'        'corr'          0.8099     0.79276     0.82576
    'reg'                '(Intercept)'        'corr'         0.86106     0.41427     0.97344
    'T'                  'T'                  'std'           2.1672      1.2363       3.799
    'siteN'              'T'                  'corr'        -0.62575    -0.66279    -0.58566
    'reg'                'T'                  'corr'        -0.89852    -0.91339    -0.88124
    'siteN'              'siteN'              'std'          0.69862      0.3096      1.5764
    'reg'                'siteN'              'corr'         0.50531    -0.29114     0.88805
    'reg'                'reg'                'std'           7.7394      3.3933      17.652

Group: Error
    Name             Estimate    Lower     Upper 
    'Res Std'        31.42       30.645    32.214
%}
