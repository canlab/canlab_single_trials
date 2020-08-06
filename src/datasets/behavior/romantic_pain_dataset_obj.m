clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2
dataset_obj = canlab_dataset;

dataset_obj.Description.Experiment_Name = 'romantic';
dataset_obj.Description.references = {'Lopez-Sola, et al. (2019) Pain'};

fprintf(['Preparing ' dataset_obj.Description.Experiment_Name '  dataset obj...\n']);

%% update these paths as needed
dataRoot     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/romantic_pain/';
% single trial image metadata
imgmeta        = importdata([dataRoot '/romantic_pain.mat']); 
% behavioral meta (might be identical to the above)
evmeta = imgmeta;

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

    n = length(evmeta.ratings{i});
    these_cov = zeros(n,length(dataset_obj.Event_Level.names));
    these_cov(:,ismember(dataset_obj.Event_Level.names,'rating')) = evmeta.ratings{i};
    these_cov(:,ismember(dataset_obj.Event_Level.names,'T')) = evmeta.temp{i};
    high_vif = zeros(n,1);
    high_vif(evmeta.high_vif_trials_idx{i}) = 1;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'high_vif')) = high_vif;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'overallN')) = 1:n;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'siteN')) = kron(ones(1,4),1:4)';
    these_cov(:,ismember(dataset_obj.Event_Level.names,'runN')) = kron(1:4,ones(1,4))';
    these_cov(:,ismember(dataset_obj.Event_Level.names,'handholding')) = [zeros(4,1);ones(8,1);zeros(4,1)];
    
    if size(these_cov,1) ~= size(img.dat,2)
        error('Metadata length mismatch');
    else
        dataset_obj.Event_Level.data{end+1} = these_cov;
    end
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
cmat = covariates;
for i = 1:length(sids)
    this_sid = sids(i);
    pain_idx = find(this_sid == covariates(:,ismember(covariate_names,'sid')));
    mean_v = mean(covariates(pain_idx,:));
    column_idx = find(ismember(covariate_names,{'rating','T','siteN','overallN'}));
    cmat(pain_idx,column_idx) = covariates(pain_idx,column_idx) - repmat(mean_v(column_idx),length(pain_idx),1);
end

tbl = array2table(cmat,'VariableNames',covariate_names);
m = fitlme(tbl,'rating ~ handholding + siteN - 1 + (handholding + siteN - 1| sid) + (handholding + siteN | runN)','FitMethod','REML')

%{

m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations             480
    Fixed effects coefficients           2
    Random effects coefficients         72
    Covariance parameters               10

Formula:
    rating ~ siteN + handholding + (siteN + handholding | sid) + (1 + siteN + handholding | runN)

Model fit statistics:
    AIC       BIC       LogLikelihood    Deviance
    3550.7    3600.7    -1763.3          3526.7  

Fixed effects coefficients (95% CIs):
    Name                 Estimate    SE         tStat      DF     pValue        Lower      Upper    
    'siteN'              -3.2021     0.53301    -6.0076    478    3.7353e-09    -4.2494      -2.1548
    'handholding'        -1.4449     0.69699    -2.0731    478      0.038697    -2.8145    -0.075397

Random effects covariance parameters (95% CIs):
Group: sid (30 Levels)
    Name1                Name2                Type          Estimate    Lower       Upper   
    'siteN'              'siteN'              'std'           1.9819      1.1602      3.3854
    'handholding'        'siteN'              'corr'        -0.70102    -0.94391    0.034091
    'handholding'        'handholding'        'std'           1.7466     0.58456      5.2188

Group: runN (4 Levels)
    Name1                Name2                Type          Estimate    Lower      Upper  
    '(Intercept)'        '(Intercept)'        'std'          2.2477     0.78303      6.452
    'siteN'              '(Intercept)'        'corr'              1         NaN        NaN
    'handholding'        '(Intercept)'        'corr'             -1         NaN        NaN
    'siteN'              'siteN'              'std'         0.26932     0.20305    0.35721
    'handholding'        'siteN'              'corr'             -1         NaN        NaN
    'handholding'        'handholding'        'std'          2.4812     0.57243     10.754

Group: Error
    Name             Estimate    Lower     Upper 
    'Res Std'        9.2753      8.6716    9.9211

%}