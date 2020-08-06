clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2
dataset_obj = canlab_dataset;

dataset_obj.Description.Experiment_Name = 'stephan';
dataset_obj.Description.references = {'Geunter, et al. (2013) Neuroimage'};

fprintf(['Preparing ' dataset_obj.Description.Experiment_Name '  dataset obj...\n']);

%% update these paths as needed
dataRoot     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/placebo_value_stephan/';
% single trial image metadata
imgmeta        = importdata([dataRoot '/placebo_value_metadata.mat']); 
% behavioral meta (might be identical to the above)
evmeta = imgmeta;

%% ideally these will match across studies
dataset_obj.Subj_Level.name = {'age','male','race','right_handed'};
dataset_obj.Subj_Level.type = {'int','boolean','level','boolean'};
dataset_obj.Event_Level.names = {'rating', 'T', 'soundintensity', 'cue', ...
    'social', 'placebo', 'value', 'runN', 'siteN', 'overallN', 'vif', ...
    'ctrl', 'reveal','conditioningN','reg','handholding','drug','open'};
dataset_obj.Event_Level.type = ...
    repmat('numeric',length(dataset_obj.Description.Event_Level));
dataset_obj.Event_Level.units = ...
    {'vas','C','levels','levels','levels','boolean','levels','count','count',...
    'count','unitless','boolean','boolean','count','levels','boolean','boolean','boolean'};

%% some study specific stuff
datset_obj.Event_Level.descrip{ismember(dataset_obj.Event_Level.names,'value')} = ...
    '-0.5 = Cheap/Weak, 0.5 = Expensive/Strong';
% see imgmeta.value_describ for source coding

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
    
    % one of the image files is missing some data, and this catches
    % corresponding entries of metadata
    this_ratings = evmeta.painrating{i};
    hasData = find(~isnan(evmeta.painrating{i}));
    n = length(hasData);
    
    these_cov = zeros(length(this_ratings(hasData)),length(dataset_obj.Event_Level.names));
    these_cov(:,ismember(dataset_obj.Event_Level.names,'rating')) = this_ratings(hasData);
    these_cov(:,ismember(dataset_obj.Event_Level.names,'T')) = imgmeta.temperature{1}(1:n); % the others are all nans
    runN = kron(1:4,ones(1,15))'; 
    these_cov(:,ismember(dataset_obj.Event_Level.names,'runN')) = runN(hasData);
    siteN = kron(ones(1,4),1:15)'; 
    these_cov(:,ismember(dataset_obj.Event_Level.names,'siteN')) = siteN(hasData);
    %overallN = kron([1,1],1:30)'; overallN = overallN(hasdata);
    % order of stimuli was resorted, so we don't have overallN info, replace with mean.
    these_cov(:,ismember(dataset_obj.Event_Level.names,'overallN'))  = 15.5*ones(n,1); 
    these_cov(:,ismember(dataset_obj.Event_Level.names,'placebo'))  = (evmeta.placebotrial{i}(hasData) + 1)/2;
    these_cov(:,ismember(dataset_obj.Event_Level.names,'value'))  = evmeta.valuetrial{i}(hasData)/2;
    %imgIdx = cellfun(@(x1)str2double(x1(6:9)),evmeta.beta_img{i});
    these_cov(:,ismember(dataset_obj.Event_Level.names,'vif')) = evmeta.vif{i}(:);
    
    if size(these_cov,1) ~= size(img.dat,2)
        error('Metadata length mismatch');
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
% center placebo
covariates(:,ismember(dataset_obj.Event_Level.names,'placebo')) = ...
    covariates(:,ismember(dataset_obj.Event_Level.names,'placebo')) - 0.5;
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
    pain_idx = find(this_sid == covariates(:,end));
    mean_v = mean(covariates(pain_idx,:));
    column_idx = find(ismember(covariate_names,{'rating','T','siteN','overallN'}));
    cmat(pain_idx,column_idx) = cmat(pain_idx,column_idx) - repmat(mean_v(column_idx),length(pain_idx),1);
end

tbl = array2table(double(cmat),'VariableNames',covariate_names);
m = fitlme(tbl,'rating ~ placebo*value + siteN - 1 + (placebo*value + siteN -1 | sid)','FitMethod','REML')
%{
m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations            2386
    Fixed effects coefficients           4
    Random effects coefficients        160
    Covariance parameters               11

Formula:
    rating ~ siteN + placebo*value + (siteN + placebo*value | sid)

Model fit statistics:
    AIC      BIC      LogLikelihood    Deviance
    23355    23442    -11663           23325   

Fixed effects coefficients (95% CIs):
    Name                   Estimate    SE        tStat      DF      pValue        Lower      Upper  
    'placebo'              -9.1938     1.8914    -4.8609    2382    1.2445e-06    -12.903    -5.4849
    'value'                -6.5434     2.1771    -3.0055    2382     0.0026788    -10.813    -2.2742
    'siteN'                 7.1867     0.2244     32.026    2382             0     6.7466     7.6267
    'placebo:value'        -7.9877     3.3662    -2.3729    2382      0.017728    -14.589    -1.3867

Random effects covariance parameters (95% CIs):
Group: sid (40 Levels)
    Name1                  Name2                  Type          Estimate    Lower       Upper  
    'placebo'              'placebo'              'std'           9.3284      6.4788     13.431
    'value'                'placebo'              'corr'        0.091951    -0.37062    0.51797
    'siteN'                'placebo'              'corr'         0.28985    -0.11288    0.61078
    'placebo:value'        'placebo'              'corr'         0.11816    -0.42335    0.59747
    'value'                'value'                'std'           11.551      8.4191     15.849
    'siteN'                'value'                'corr'        -0.18091    -0.50752    0.19116
    'placebo:value'        'value'                'corr'         0.68914    0.069508    0.92506
    'siteN'                'siteN'                'std'           1.4192      1.1367     1.7718
    'placebo:value'        'siteN'                'corr'        -0.25282    -0.61174    0.19244
    'placebo:value'        'placebo:value'        'std'           15.131      9.7362     23.516

Group: Error
    Name             Estimate    Lower     Upper
    'Res Std'        28.889      28.053    29.75
%}
