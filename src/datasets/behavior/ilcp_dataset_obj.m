% ILCP could use some more nuanced updates to the cue paradigm. Right now
% it indicates what the subject expects, which is a composite of some cues
% and some probabilities. It could be sutdivided further to indicate the
% actual cue (low or high) and the actual probabilities (80/20 vs 5/50) 
% which are presented to the subject.

clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2

fprintf('Preparing ILCP dataset obj...\n');
dataset_obj = canlab_dataset;
dataset_obj.Description.Experiment_Name = 'ilcp';
dataset_obj.Description.references = {'Woo, et al. (2017) Nature Communications'};

%% update these paths as needed
dataRoot     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/ILCP_wani/';
% single trial image metadata
imgmeta        = importdata([dataRoot '/ilcp_metainfo_waniupdate.mat']); 
% behavioral meta (might be identical to the above)
evmeta = importdata('/work/ics/data/projects/wagerlab/labdata/current/bogdan_spillover_storage/Wagerlab_Single_Trial_Pain_Datasets/Data/ILCP_wani/ilcp_metainfo_waniupdate.mat');

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
        
    ref_ratings = imgmeta.ratings{i};
    meta_idx = find(~isnan(evmeta.ratings{i}));
    if length(imgmeta.ratings{i}) ~= length(evmeta.ratings{i}) || length(imgmeta.ratings{i}) ~= size(img.dat,2)
        error('Metadata length mismatch');
    end

    % looks like behavioral covariates are contained here:
    % /work/ics/data/projects/wagerlab/labdata/current/ilcp/Behavioral/Sub29/PainTestSub29Session3Run2.mat
    % there are several manipulations, has/lacks control and low/high
    % expectation cues, but then one or the other cue can be selected, so
    % we have 8 levels, and one of them should be under represented (has
    % control, high expectation conditions are unlikely to lead to the
    % selection of the cue with high probability of high pain).
    % The data at the above path seems to have 4 columns though, instead of
    % 2, and I don't understand what the difference is between "choice"
    % (presumably choice of cue) and "response". Nor can I tell which
    % numbers correspond to which of the two conditions (not that I
    % necessarily need to be able to tell for coding my models, but does 
    % help for sanity checking). 

    fnum_idx = strfind(evmeta.subjects{i},'_') - 1;
    fnum = evmeta.subjects{i}(5:fnum_idx(1));

    subj_path = ['/work/ics/data/projects/wagerlab/labdata/current/ilcp/Behavioral/*ub', fnum, '/'];
    choice_fnames = dir([subj_path, '/PainTestSub', fnum, 'Session*Run*.mat']);
    obs_fnames = dir([subj_path, '/obsPainTestSub', fnum, 'Session*Run*.mat']);

    subj_dat = [];
    for k = 1:length(choice_fnames)
        this_dat = importdata([choice_fnames(k).folder, '/', choice_fnames(k).name]);
        this_dat = [this_dat.paindata, ones(size(this_dat.paindata,1),1)];
        subj_dat = [subj_dat; this_dat];
    end
    for k = 1:length(obs_fnames)
        this_dat = importdata([obs_fnames(k).folder, '/', obs_fnames(k).name]);
        this_dat = [this_dat.obspaindata(:,1:end-1), zeros(size(this_dat.obspaindata,1),1)];
        subj_dat = [subj_dat; this_dat];
    end

    subj_dat = [subj_dat(:,1),kron(1:8,ones(1,8))',subj_dat(:,2:end)];
    subj_dat = sortrows(subj_dat,1); % this aligns the data with the single trial data sequence
    subj_dat = subj_dat(:,2:end);
    % for labels see
    % /work/ics/data/projects/wagerlab/labdata/current/ilcp/Behavioral/GetBehavDataPR.m
    % Note that "response" is left or right mouse button. Choice indicates
    % the "cue". See Liane's email from Nov 21st 2018
    subj_dat = array2table(subj_dat,'VariableNames',{'run','trial','pair','choice','response','control','expected_pain','felt_pain','has_control'});
    subj_dat = subj_dat(~isnan(subj_dat.felt_pain),:);

    if corr(subj_dat.felt_pain,ref_ratings) + 0.0001 < 1 || corr(evmeta.ratings{i}(meta_idx),ref_ratings) + 0.0001 < 1
        error(['Could not import behavioral data for subject ' num2str(i)']);
    end

    ratings = evmeta.ratings{i}(meta_idx);
    temp = evmeta.temp{i}(meta_idx);

    % The coding scheme was deduced from 
    % /work/ics/data/projects/wagerlab/labdata/current/ilcp/Behavioral/GetBehavData.PR.m
    cue = zeros(length(ratings),1);
    idx_80_20 = find(subj_dat.pair == 1); % 80/20 cue pair
    idx_50 = find(subj_dat.pair == 2);

    % double check that these should be "choice" and not "response"
    idx_20 = idx_80_20(subj_dat.choice(idx_80_20) == -1); % 20 cue
    idx_80 = idx_80_20(subj_dat.choice(idx_80_20) == 1); % 80 cue
    idx_50a = idx_50(subj_dat.choice(idx_50) == -1);
    idx_50b = idx_50(subj_dat.choice(idx_50) == 1);

    cue(idx_20) = 3;
    cue(idx_50a) = 2.5;
    cue(idx_50b) = 1.5;
    cue(idx_80) = 1;
    cue(cue == 0) = nan;

    ctrl = subj_dat.has_control;
    runN = subj_dat.run;
    siteN = subj_dat.trial;
    overallN = meta_idx;
    vif = double(evmeta.high_vif_trials_idx{i}(meta_idx));
    vif(vif > 0) = 4;

    these_cov = zeros(length(ratings),length(dataset_obj.Event_Level.names));
    these_cov(:,contains(dataset_obj.Event_Level.names,'rating')) = ratings;
    these_cov(:,contains(dataset_obj.Event_Level.names,'T')) = temp;
    these_cov(:,contains(dataset_obj.Event_Level.names,'cue')) = cue - 2;
    these_cov(:,contains(dataset_obj.Event_Level.names,'ctrl')) = ctrl;
    these_cov(:,contains(dataset_obj.Event_Level.names,'runN')) = runN;
    these_cov(:,contains(dataset_obj.Event_Level.names,'siteN')) = siteN;
    these_cov(:,contains(dataset_obj.Event_Level.names,'overallN')) = overallN;
    these_cov(:,contains(dataset_obj.Event_Level.names,'vif')) = vif;
    
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
uniq_sid = unique(sids);
cmat = covariates;
for j = 1:length(uniq_sid)
    this_sid = uniq_sid(j);
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
tbl.cue = tbl.cue - 2;
tbl.ctrl = tbl.ctrl - 0.5;
m = fitlme(tbl,'rating ~ T + cue*ctrl - 1 + (T + cue*ctrl - 1| sid) + (T | runN)','FitMethod','REML')

%{
m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations            1826
    Fixed effects coefficients           4
    Random effects coefficients        132
    Covariance parameters               14

Formula:
    rating ~ T + cue*ctrl + (T + cue*ctrl | sid) + (1 + T | runN)

Model fit statistics:
    AIC      BIC      LogLikelihood    Deviance
    14414    14514    -7189.2          14378   

Fixed effects coefficients (95% CIs):
    Name              Estimate    SE         tStat       DF      pValue        Lower      Upper  
    'T'                 11.274    0.92519      12.185    1822             0     9.4593     13.088
    'cue'              0.48659    0.39576      1.2295    1822       0.21904    -0.2896     1.2628
    'ctrl'             -3.6723     1.0102     -3.6353    1822    0.00028537    -5.6535    -1.6911
    'cue:ctrl'        -0.32412    0.81589    -0.39727    1822       0.69122    -1.9243      1.276

Random effects covariance parameters (95% CIs):
Group: sid (29 Levels)
    Name1             Name2             Type          Estimate    Lower       Upper   
    'T'               'T'               'std'           4.0268      3.0094      5.3881
    'cue'             'T'               'corr'        -0.33649         NaN         NaN
    'ctrl'            'T'               'corr'        -0.61942    -0.80263    -0.32941
    'cue:ctrl'        'T'               'corr'         0.79667         NaN         NaN
    'cue'             'cue'             'std'          0.85678     0.40834      1.7977
    'ctrl'            'cue'             'corr'        -0.12055         NaN         NaN
    'cue:ctrl'        'cue'             'corr'        -0.19954         NaN         NaN
    'ctrl'            'ctrl'            'std'           4.4794       3.063      6.5507
    'cue:ctrl'        'ctrl'            'corr'        -0.94074    -0.97246     -0.8748
    'cue:ctrl'        'cue:ctrl'        'std'           2.0211     0.88365      4.6228

Group: runN (8 Levels)
    Name1                Name2                Type          Estimate    Lower       Upper 
    '(Intercept)'        '(Intercept)'        'std'         0.23133     0.003113     17.19
    'T'                  '(Intercept)'        'corr'              1          NaN       NaN
    'T'                  'T'                  'std'          1.3679      0.67832    2.7586

Group: Error
    Name             Estimate    Lower     Upper 
    'Res Std'        12.091      11.698    12.497

%}

tbl_raw = array2table(double(covariates),'VariableNames',covariate_names);
fprintf('The below should reproduce figure 5c in the SIIPS paper:\n');
HC = tbl_raw.ctrl == 1;
LC = tbl_raw.ctrl == 0;
LE = (tbl_raw.cue == 1 | tbl_raw.cue == 3); %80/20 pairs, presumably subject is smart and picks 80 mostly, hence "low" expectency"
HE = tbl_raw.cue == 1.5 | tbl_raw.cue == 2.5; %50/50 pairs

[nanmean(tbl_raw.rating((HC .* LE) == 1)), ...
nanmean(tbl_raw.rating((HC .* HE) == 1)), ...
nanmean(tbl_raw.rating((LC .* LE) == 1)), ...
nanmean(tbl_raw.rating((LC .* HE) == 1))]

%{
ans =

   32.6386   34.7808   35.3470   38.9526
%}