clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2

fprintf('Preparing BMRK4 dataset obj...\n');
dataset_obj = canlab_dataset;
dataset_obj.Description.Experiment_Name = 'bmrk4';
dataset_obj.Description.references = {'Krishnan et al. (2016) ELife'};

%% update these paths as neededf
dataRoot = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/bmrk4_smoothed_withbasis/';
% single trial image metadata
imgmeta        = importdata([dataRoot '/bmrk4_single_trial_model.mat']); 
% behavioral meta (might be identical to the above)
evmeta = importdata('/work/ics/data/projects/wagerlab/labdata/single_trials_experimental/Data/bmrk4_smoothed_withbasis/bmrk4_single_trial_model.mat');
evDataRoot = '/work/ics/data/projects/wagerlab/labdata/current/bmrk4/Data/BMRK_Data/';

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

    % check for mismatches in ev and image data. 
    evnum = ~isnan(evmeta.ratings{i});
    imgnum = ~isnan(imgmeta.ratings{i});
    if corr(evmeta.ratings{i}(evnum),imgmeta.ratings{i}(imgnum)) + 0.001 < 1
        %warning(['Skipping subject ' int2str(i) ' due to data length mismatch']);
        %continue;
        error(['subject ' int2str(i) ' length mismatch']);
    end
    
    % import covariates from evmeta file(s)
    this_subj = evmeta.subjects{i};
    numstr = this_subj(4:end);
    subjRoot = [evDataRoot, 'sub_', numstr, '/Pain/'];
    fnames = dir([subjRoot, 'bmrk4_env_sub', numstr, '_Pain_run*mat']);
    these_cov = [];
    % iterate through one run at a time
    for k = 3:length(fnames) %skip first two, these are training runs and not included in single trial dataset
        M = importdata([subjRoot, fnames(k).name]);
        loc = M.info.painlocation;
        seq1 = 1:length(M.rating);
        next_cov = zeros(length(M.rating),length(dataset_obj.Event_Level.names));
        next_cov(:,ismember(dataset_obj.Event_Level.names,...
            {'rating','T','cue','runN','siteN'})) = ...
        	[M.rating, M.stim, M.cue - 2, repmat(loc,length(M.rating),1),seq1(:)];
        idx = cell(length(M.TEMPERATURES),1);
        % swap temperature levels for actual intensities
        for l = 1:length(M.TEMPERATURES)
            idx{l} = find(next_cov(:,ismember(dataset_obj.Event_Level.names,'T')) == l);
        end
        for l = 1:length(M.TEMPERATURES)
            next_cov(idx{l},ismember(dataset_obj.Event_Level.names,'T')) = M.TEMPERATURES(l);
        end
        these_cov = [these_cov; next_cov];
    end
    seq2 = 1:size(these_cov,1);
    vif = zeros(size(these_cov,1),1);
    vif(evmeta.high_vif_trials_idx{i}) = 1;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'overallN')) = seq2;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'high_vif')) = vif;
    
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
sids = unique(sids);
cmat = covariates;
for i = 1:length(sids)
    this_sid = sids(i);
    pain_idx = find(this_sid == covariates(:,ismember(covariate_names,'sid')));
    mean_v = mean(covariates(pain_idx,:));
    column_idx = find(ismember(covariate_names,{'rating','T','siteN','overallN'}));
    cmat(pain_idx,column_idx) = cmat(pain_idx,column_idx) - repmat(mean_v(column_idx),length(pain_idx),1);
end
cmat(:,contains(covariate_names,'cue')) = cmat(:,contains(covariate_names,'cue'))-2;

covariates(:,contains(covariate_names,'cue')) = covariates(:,contains(covariate_names,'cue')) - 2;

%try
%    dat{1}.covariates = covariates;
%    dat{1}.covariate_names = covariate_names;
%end

tbl = array2table(double(cmat),'VariableNames',covariate_names);
m = fitlme(tbl,'rating ~ T + cue + siteN - 1 + (T + cue + siteN - 1| sid) + (T + cue + siteN | runN)','FitMethod','REML')

%{

m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations            2268
    Fixed effects coefficients           3
    Random effects coefficients        116
    Covariance parameters               17

Formula:
    rating ~ T + cue + siteN + (T + cue + siteN | sid) + (1 + T + cue + siteN | runN)

Model fit statistics:
    AIC      BIC      LogLikelihood    Deviance
    16478    16592    -8218.9          16438   

Fixed effects coefficients (95% CIs):
    Name           Estimate    SE         tStat      DF      pValue        Lower       Upper    
    'T'              10.919    0.88536     12.333    2265             0      9.1825       12.655
    'cue'            2.4099    0.32825     7.3416    2265    2.9288e-13      1.7662       3.0536
    'siteN'        -0.21317    0.10123    -2.1059    2265      0.035324    -0.41168    -0.014666

Random effects covariance parameters (95% CIs):
Group: sid (28 Levels)
    Name1          Name2          Type          Estimate    Lower       Upper  
    'T'            'T'            'std'           4.4415      3.3357     5.9138
    'cue'          'T'            'corr'         0.39132    -0.17756    0.76418
    'siteN'        'T'            'corr'        -0.23221    -0.67035      0.326
    'cue'          'cue'          'std'           1.2041     0.70687     2.0512
    'siteN'        'cue'          'corr'        -0.21855     -0.7648    0.51052
    'siteN'        'siteN'        'std'          0.36897     0.21276    0.63985

Group: runN (8 Levels)
    Name1                Name2                Type          Estimate    Lower         Upper   
    '(Intercept)'        '(Intercept)'        'std'           1.2114       0.68281      2.1493
    'T'                  '(Intercept)'        'corr'         0.44685       0.43963     0.45401
    'cue'                '(Intercept)'        'corr'        -0.33839           NaN         NaN
    'siteN'              '(Intercept)'        'corr'         0.93368       0.93036     0.93684
    'T'                  'T'                  'std'          0.51436        0.1502      1.7615
    'cue'                'T'                  'corr'        -0.99304           NaN         NaN
    'siteN'              'T'                  'corr'        0.096848           NaN         NaN
    'cue'                'cue'                'std'          0.19878      0.010133      3.8995
    'siteN'              'cue'                'corr'        0.021034    -0.0029565    0.045001
    'siteN'              'siteN'              'std'         0.083246      0.012711     0.54521

Group: Error
    Name             Estimate    Lower     Upper 
    'Res Std'        8.8182      8.5598    9.0844

%}