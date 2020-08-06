clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2

fprintf('Preparing IE2 dataset obj...\n');
dataset_obj = canlab_dataset;
dataset_obj.Description.Experiment_Name = 'ie2';
dataset_obj.Description.references = {'Jepma et al. (2018) Nature Human Behavior'};


%% update these paths as needed
dataRoot     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/IE2_NEW/';
% single trial image metadata
imgmeta        = importdata([dataRoot '/ie2_metadata.mat']); 
% behavioral meta (might be identical to the above)
evmeta1 = importdata('/work/ics/data/projects/wagerlab/labdata/single_trials_experimental/Data/IE2/ie2_metadata.mat');
evmeta2 = importdata('/work/ics/data/projects/wagerlab/labdata/current/ie2/Expect_CS_Pain_Trial_Temp_PE_Sub.mat');

[demoNum, demoTxt] = xlsread('demographics_IE2_fmri.xlsx');


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
dataset_obj.Subj_Level.descrip{ismember(dataset_obj.Subj_Level.names,'race')} = ...
    '1: white (not hispanic origin), 2: hispanic, 3: Asian or pacific islander';
datset_obj.Event_Level.descrip{ismember(dataset_obj.Event_Level.names,'cue')} = ...
    ['Cue actually has 7 levels, the positive and negative ones (3 ea.) are ',...
    'collapsed here. See ie2_datset_obj.m (dataset_obj generating script)', ...
    'for details. Theres a linear model commented out there which',...
    'illustrates the variable effects across the 7 levels.'];
datset_obj.Event_Level.descrip{ismember(dataset_obj.Event_Level.names,'high_vif')} = ...
    'vif > 3.0';

% The following code  evaluates the effect sizes of each cue to identify
% which are low cues, which are high cues, and which is the neutral cue
% If you run it you'll find cue 1-3 are low, 4-6 are high and 7 is neutral.
% They have different effect sizes though (so the 3 low are distinguishable
% from one another, ditto for the 3 high)
%{
x = cell2mat(evmeta2');

cue = [];
for i=1:7
    cue(:,i) = (x(:,2) == i);
end

design = [x(:,3), x(:,5), cue(:,1:7), x(:,end)];
sid = unique(design(:,end));
for i = 1:length(sid) %subject center design
    this_sid = sid(i);
    idx = find(this_sid == design(:,end));
    design(idx,1:2) = design(idx,1:2) - repmat(mean(design(idx,1:2)),length(idx),1);
end
tbl = array2table(design,'VariableNames',{'pain','T','cue1','cue2','cue3','cue4','cue5','cue6','cue7','sid'});

m = fitlme(tbl,'pain ~ T + cue1 + cue2 + cue3 + cue4 + cue5 + cue6 + cue7 - 1 + (T + cue1 + cue2 + cue3 + cue4 + cue5 + cue6 + cue7 - 1| sid)','FitMethod', 'REML');
%}

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

    ev1_idx = find(~isnan(evmeta1.ratings{i}));
    cue_mat_ratings = evmeta2{i}(:,3);
    siteN = kron(ones(1,5),1:14)'; 
    runN = kron(1:5,ones(1,14))';
    
    ref_ratings = imgmeta.ratings{i};
    img_idx = find(~isnan(ref_ratings));

    these_cov = [];
    if corr(ref_ratings(img_idx), cue_mat_ratings) + 0.00001 > 1 && length(evmeta1.ratings{i}) == size(img.dat,2)
        these_cov = zeros(length(ref_ratings),length(dataset_obj.Event_Level.names));
        these_cov(:,ismember(dataset_obj.Event_Level.names,'rating')) = evmeta2{i}(:,3);
        these_cov(:,ismember(dataset_obj.Event_Level.names,'T')) = evmeta2{i}(:,5);
        these_cov(:,ismember(dataset_obj.Event_Level.names,'cue')) = evmeta2{i}(:,2);
        these_cov(:,ismember(dataset_obj.Event_Level.names,'runN')) = runN;
        these_cov(:,ismember(dataset_obj.Event_Level.names,'siteN')) = siteN;
        these_cov(:,ismember(dataset_obj.Event_Level.names,'overallN')) = evmeta2{i}(:,4);

        % recode labels according to estimated cue identity. For cue
        % identity estimation see commented LME code above.
        cue = cell(7,1);
        for k = 1:7
            cue{k} = find(these_cov(:,ismember(dataset_obj.Event_Level.names,'cue')) == k);
        end
        for k = 1:3
            these_cov(cue{k},ismember(dataset_obj.Event_Level.names,'cue')) = -0.5;
        end
        for k = 4:6
            these_cov(cue{k},ismember(dataset_obj.Event_Level.names,'cue')) = 0.5;
        end
        these_cov(cue{7},ismember(dataset_obj.Event_Level.names,'cue')) = 0;
    end

    nonnan = ~isnan(evmeta1.ratings{i}); 
    if corr(ref_ratings(img_idx),evmeta1.ratings{i}(nonnan)) + 0.00001 > 1
        these_cov(evmeta1.high_vif_trials_idx{i} == 1, contains(dataset_obj.Event_Level.names,'high_vif')) = 1;
    else
        warning(['Could not import VIF values for subject ' num2str(i)]);
    end
    
    dataset_obj.Event_Level.data{end+1} = these_cov;
    
    % get demographic info
    these_demos = zeros(1,length(dataset_obj.Subj_Level.names));
    ageIdx = find(contains(demoTxt(1,:)','Age'));
    raceIdx = find(contains(demoTxt(1,:)','Race'));
    genderIdx = find(contains(demoTxt(1,:)','Gender'));
    handIdx = find(contains(demoTxt(1,:)','Hand'));
    
    demoIdx = find(arrayfun(@(x1)(contains(subj,num2str(x1))),demoNum(:,3)));
    if isempty(demoIdx)
        warning(['Could not retrieve matching demographic info for subject ' int2str(i)]);
        these_demos(:) = nan;
    else

        if ~isempty(demoNum(demoIdx,ageIdx)) && ~isnan(demoNum(demoIdx,ageIdx))
            these_demos(ismember(dataset_obj.Subj_Level.names,'age')) = demoNum(demoIdx,ageIdx);
        else
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
uniq_sid = unique(sids);
cmat = covariates;
for i = 1:length(uniq_sid)
    this_sid = uniq_sid(i);
    pain_idx = find(this_sid == covariates(:,ismember(covariate_names,'sid')));
    mean_v = mean(covariates(pain_idx,:));
    column_idx = find(ismember(covariate_names,{'rating','T','siteN','overallN'}));
    cmat(pain_idx,column_idx) = cmat(pain_idx,column_idx) - repmat(mean_v(column_idx),length(pain_idx),1);
end

tbl = array2table(cmat,'VariableNames',covariate_names);
m = fitlme(tbl,'rating ~ T + cue + siteN - 1 + (T + cue + siteN - 1| sid) + (T + siteN | runN)','FitMethod','REML')
%{
m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations            1330
    Fixed effects coefficients           3
    Random effects coefficients         72
    Covariance parameters               13

Formula:
    rating ~ T + cue + siteN + (T + cue + siteN | sid) + (1 + T + siteN | runN)

Model fit statistics:
    AIC       BIC     LogLikelihood    Deviance
    9388.9    9472    -4678.5          9356.9  

Fixed effects coefficients (95% CIs):
    Name           Estimate    SE          tStat      DF      pValue        Lower       Upper   
    'T'              4.2421     0.67923     6.2454    1327    5.6801e-10      2.9096      5.5746
    'cue'             9.041      1.3992     6.4617    1327    1.4501e-10      6.2962      11.786
    'siteN'        -0.13846    0.089722    -1.5432    1327       0.12303    -0.31447    0.037556

Random effects covariance parameters (95% CIs):
Group: sid (19 Levels)
    Name1          Name2          Type          Estimate     Lower       Upper  
    'T'            'T'            'std'            2.2642      1.2958     3.9565
    'cue'          'T'            'corr'        -0.039448     -0.5916     0.5379
    'siteN'        'T'            'corr'         -0.15455    -0.56644    0.31914
    'cue'          'cue'          'std'            5.7416      3.9713      8.301
    'siteN'        'cue'          'corr'         -0.08439    -0.66913    0.56489
    'siteN'        'siteN'        'std'           0.21784     0.11068    0.42875

Group: runN (5 Levels)
    Name1                Name2                Type          Estimate    Lower         Upper  
    '(Intercept)'        '(Intercept)'        'std'         0.67023        0.29395     1.5282
    'T'                  '(Intercept)'        'corr'              1            NaN        NaN
    'siteN'              '(Intercept)'        'corr'             -1            NaN        NaN
    'T'                  'T'                  'std'         0.11899     5.7302e-05     247.08
    'siteN'              'T'                  'corr'             -1            NaN        NaN
    'siteN'              'siteN'              'std'         0.19363       0.082895    0.45231

Group: Error
    Name             Estimate    Lower     Upper 
    'Res Std'        7.9307      7.6282    8.2452
%}