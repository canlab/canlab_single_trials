clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2
dataset_obj = canlab_dataset;

dataset_obj.Description.Experiment_Name = 'ie';
dataset_obj.Description.references = {'Roy, et al. (2014) Nature Neuroscience'};

fprintf(['Preparing ' dataset_obj.Description.Experiment_Name '  dataset obj...\n']);

%% update these paths as needed
dataRoot     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/ie_for_tor/';
% single trial image metadata
imgmeta        = importdata([dataRoot '/ie_model.mat']); 
% behavioral meta (might be identical to the above)
evmeta1 = importdata('/work/ics/data/projects/wagerlab/labdata/single_trials_experimental/Data/ie_for_tor/ie_model.mat');
evmeta2 = importdata('/projects/bope9760/pain_predictive_rois_rnd4/resources/behavioral_analyses/ie_metadata.mat');

[demoNum, demoTxt] = xlsread('IE_Demographics_fMRI.xlsx');

%% ideally these will match across studies
dataset_obj.Subj_Level.names = {'age','gender','race','right_handed'};
dataset_obj.Subj_Level.type = {'int','level','level','boolean'};
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

    imgidx = ~isnan(imgmeta.ratings{i});
    evidx = ~isnan(evmeta1.ratings{i});
    n = length(imgmeta.ratings{i});
    if n == size(img.dat,2) && corr(imgmeta.ratings{i}(imgidx),evmeta1.ratings{i}(evidx))+0.001 > 1
        these_cov = zeros(n,length(dataset_obj.Event_Level.names ));
        these_cov(:,ismember(dataset_obj.Event_Level.names ,'rating')) = imgmeta.ratings{i};
        these_cov(:,ismember(dataset_obj.Event_Level.names ,'T')) = imgmeta.temp{i} + 47;
        these_cov(:,ismember(dataset_obj.Event_Level.names ,{'overallN','siteN'})) = mean(1:n)*ones(n,2);
        these_cov(:,ismember(dataset_obj.Event_Level.names ,{'runN'})) = ones(n,1);
        cues = evmeta1.cues{i};
        cues(cues == 1 | cues == 3) = -0.5;
        cues(cues == 2 | cues == 4) = 0.5;
        these_cov(:,ismember(dataset_obj.Event_Level.names ,'cue')) = cues;
        these_cov(:,ismember(dataset_obj.Event_Level.names ,'placebo')) = evmeta2.placebo{i}/2 + 0.5;
    else
        error('Metadata length mismatch');
    end

    
    dataset_obj.Event_Level.data{end+1} = these_cov;
    
    % get subject level data
    
    these_demos = zeros(1,length(dataset_obj.Subj_Level.names));
    % find age based on column that has
    % a) all entries < 100
    % b) all entries positive
    % c) all entries are round numbers
    ageIdx = find(all(demoNum < 100) & all(demoNum > 0) & all(arrayfun(@(x1)(x1 == int32(x1)),demoNum)));
    raceIdx = find(contains(demoTxt(1,:)','race'));
    genderIdx = find(contains(demoTxt(1,:)','gender'));
    handIdx = find(contains(demoTxt(1,:)','Handedeness?'));
    
    demoIdx = find(arrayfun(@(x1)(contains(subj,num2str(x1))),demoNum(:,8)));
    if isempty(demoIdx)
        demoIdx = find(arrayfun(@(x1)(contains(subj,num2str(x1))),demoNum(:,2)));
        if isempty(demoIdx)
            warning(['Could not retrieve matching demographic info for subject ' int2str(i)]);
            these_demos(:) = nan;
        end
    end
    if ~isempty(demoIdx)
        if ~isempty(demoNum(demoIdx,ageIdx)) && ~isnan(demoNum(demoIdx,ageIdx))
            these_demos(ismember(dataset_obj.Subj_Level.names,'age')) = demoNum(demoIdx,ageIdx);
        else
            warning(['Could not retrieve age for subject ' int2str(i) ]);
            these_demos(ismember(dataset_obj.Subj_Level.names,'age')) = nan;
        end
        
        if contains(demoTxt{demoIdx + 1,genderIdx},'Male')
            these_demos(ismember(dataset_obj.Subj_Level.names,'male')) = 1;
        elseif contains(demoTxt{demoIdx + 1,genderIdx},'Female')
            these_demos(ismember(dataset_obj.Subj_Level.names,'male')) = 0;
        else
            these_demos(ismember(dataset_obj.Subj_Level.names,'male')) = nan;
            warning(['Could not retrieve gender for subject ' int2str(i)]);
        end

        switch demoTxt{demoIdx + 1,raceIdx}
            case 'White (not of Hispanic origin)'
                these_demos(ismember(dataset_obj.Subj_Level.names,'race')) = 1;
            case 'Hispanic'
                these_demos(ismember(dataset_obj.Subj_Level.names,'race')) = 2;
            case 'Asian or Pacific Islander'
                these_demos(ismember(dataset_obj.Subj_Level.names,'race')) = 3;
            otherwise
                these_demos(ismember(dataset_obj.Subj_Level.names,'race')) = nan;
                warning(['Could not retrieve race for subject ' int2str(i) ]);
        end

        switch strrep(demoTxt{demoIdx + 1, handIdx},';','')
            case 'Right'
                these_demos(ismember(dataset_obj.Subj_Level.names,'right_handed')) = 1;
            case 'Left'
                these_demos(ismember(dataset_obj.Subj_Level.names,'right_handed')) = 0;
            otherwise
                these_demos(ismember(dataset_obj.Subj_Level.names,'right_handed')) = nan;
                warning(['Could not retrieve hand for subject ' int2str(i) ]);
        end
    end
    dataset_obj.Subj_Level.data(end+1,:) = these_demos;
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
covariates = covariates(~any(isnan(covariates),2),:);
uniq_sid = unique(covariates(:,ismember(covariate_names,'sid')));
cmat = covariates;
for i = 1:length(uniq_sid)
    this_sid = uniq_sid(i);
    pain_idx = find(this_sid == covariates(:,ismember(covariate_names,'sid')));
    mean_v = mean(covariates(pain_idx,:));
    column_idx = find(ismember(covariate_names,{'rating','T','cue','placebo','siteN','overallN'}));
    cmat(pain_idx,column_idx) = cmat(pain_idx,column_idx) - repmat(mean_v(column_idx),length(pain_idx),1);
end

tbl = array2table(cmat,'VariableNames',covariate_names);
m = fitlme(tbl,'rating ~ T + cue + placebo - 1 + (T + cue + placebo | sid)','FitMethod','REML')
%{
m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations            2400
    Fixed effects coefficients           3
    Random effects coefficients        200
    Covariance parameters               11

Formula:
    rating ~ T + cue + placebo + (1 + T + cue + placebo | sid)

Model fit statistics:
    AIC      BIC      LogLikelihood    Deviance
    17475    17556    -8723.4          17447   

Fixed effects coefficients (95% CIs):
    Name             Estimate    SE         tStat      DF      pValue        Lower      Upper  
    'T'               9.8813      0.5957     16.588    2397             0     8.7131     11.049
    'cue'             2.7758     0.55409     5.0098    2397    5.8473e-07     1.6893     3.8624
    'placebo'        -8.7779       1.289    -6.8099    2397    1.2306e-11    -11.306    -6.2503

Random effects covariance parameters (95% CIs):
Group: sid (50 Levels)
    Name1                Name2                Type          Estimate      Lower       Upper  
    '(Intercept)'        '(Intercept)'        'std'         5.7342e-11         NaN        NaN
    'T'                  '(Intercept)'        'corr'        7.9787e-07         NaN        NaN
    'cue'                '(Intercept)'        'corr'        3.0273e-05         NaN        NaN
    'placebo'            '(Intercept)'        'corr'        6.2975e-06         NaN        NaN
    'T'                  'T'                  'std'             3.3396      2.4878      4.483
    'cue'                'T'                  'corr'           0.59658    -0.58598    0.96721
    'placebo'            'T'                  'corr'          -0.61826    -0.77903    -0.3812
    'cue'                'cue'                'std'             1.4733     0.41655     5.2109
    'placebo'            'cue'                'corr'          -0.51471    -0.92719    0.46193
    'placebo'            'placebo'            'std'              4.875      2.5669     9.2586

Group: Error
    Name             Estimate    Lower    Upper 
    'Res Std'        8.8927      8.638    9.1549
%}