% note there are emotional faces that subjects are cued with, but that data
% is missing
clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2

fprintf('Preparing NSF dataset obj...\n');
dataset_obj = canlab_dataset;
dataset_obj.Description.Experiment_Name = 'nsf';
dataset_obj.Description.references = {'Wager, et al. (2013) New England Journal of Medicine', 'Atlas, et al. (2014) Pain'};

%% update these paths as needed
dataRoot     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/nsf/';
% single trial image metadata
imgmeta        = importdata([dataRoot '/meta.mat']); 
% behavioral meta (might be identical to the above)
evmeta = importdata('/work/ics/data/projects/wagerlab/labdata/single_trials_experimental/Data/NSF/meta.mat'); 

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

    these_cov = zeros(48,length(dataset_obj.Event_Level.names));
    evnum = ~isnan(evmeta.ratings{i});
    imgnum = ~isnan(imgmeta.ratings{i});
    if corr(evmeta.ratings{i}(evnum),imgmeta.ratings{i}(imgnum)) + 0.001 < 1
        warning(['Skipping subject ' int2str(i) ' due to data length mismatch']);
        continue;
    end
    indices = zeros(48,1);
    for k = 1:length(evmeta.ratings{i})
        if ~isnan(evmeta.ratings{i}(k))
            imgNumIdx = strfind(evmeta.images{i}{k},'.img');
            indices(str2num(evmeta.images{i}{k}(imgNumIdx-2:imgNumIdx-1))) = 1;
        end
    end
    these_cov(find(indices),ismember(dataset_obj.Event_Level.names,'rating')) = evmeta.ratings{i}(evnum);  % rating
    these_cov(find(indices),ismember(dataset_obj.Event_Level.names,'T')) = evmeta.temp{i}(evnum);     % T
    these_cov(:,ismember(dataset_obj.Event_Level.names,'runN')) = reshape(kron(1:6,ones(8,1)),48,1);     % runN
    these_cov(:,ismember(dataset_obj.Event_Level.names,'siteN')) = repmat(1:8,1,6)';                      % siteN
    these_cov(:,ismember(dataset_obj.Event_Level.names,'overallN')) = 1:48;                                  % overallN
    
    dataset_obj.Event_Level.data{end+1} = these_cov(indices == 1,:);
end

% drop empty fields
hasData = any(cell2mat(dataset_obj.Event_Level.data'));
for i = 1:length(dataset_obj.Event_Level.data)
    dataset_obj.Event_Level.data{i} = dataset_obj.Event_Level.data{i}(:,hasData);
end
dataset_obj.Event_Level.names = dataset_obj.Event_Level.names(hasData);
dataset_obj.Event_Level.type = dataset_obj.Event_Level.type(hasData);
dataset_obj.Event_Level.units = dataset_obj.Event_Level.units(hasData);

save([analysisRoot 'resources/prep_canlab_dataset_objs/' ...
    dataset_obj.Description.Experiment_Name, '_dataset_obj.mat'],'dataset_obj');


%% sanity check
covariates = cell2mat(dataset_obj.Event_Level.data');
% add subj id column
sids = [];
for i = 1:length(dataset_obj.Event_Level.data)
    sids = [sids; i*ones(size(dataset_obj.Event_Level.data{i},1),1)];
end
covariates = [covariates, sids];
covariate_names = [dataset_obj.Event_Level.names, 'sid'];

% center within subject
cmat = covariates;
for i = 1:length(sids)
    this_sid = sids(i);
    pain_idx = find(this_sid == covariates(:,ismember(covariate_names,'sid')));
    mean_v = mean(covariates(pain_idx,:));
    column_idx = find(ismember(covariate_names,{'rating','T','siteN','overallN'}));
    cmat(pain_idx,column_idx) = covariates(pain_idx,column_idx) - repmat(mean_v(column_idx),length(pain_idx),1);
end

tbl = array2table(cmat,'VariableNames',covariate_names);
m = fitlme(tbl,'rating ~ T + siteN - 1 + (T + siteN - 1| sid) + (T + siteN | runN)','FitMethod','REML')

%{

m = 


Linear mixed-effects model fit by REML

Model information:
    Number of observations            1149
    Fixed effects coefficients           3
    Random effects coefficients         96
    Covariance parameters               13

Formula:
    rating ~ 1 + T + siteN + (1 + T + siteN | sid) + (1 + T + siteN | runN)

Model fit statistics:
    AIC       BIC       LogLikelihood    Deviance
    4163.6    4244.3    -2065.8          4131.6  

Fixed effects coefficients (95% CIs):
    Name                 Estimate    SE          tStat      DF      pValue        Lower       Upper    
    '(Intercept)'         -28.165      2.3089    -12.198    1146             0     -32.695      -23.634
    'T'                   0.72653    0.051568     14.089    1146             0     0.62535      0.82771
    'siteN'              -0.10727    0.029949    -3.5817    1146    0.00035573    -0.16603    -0.048506

Random effects covariance parameters (95% CIs):
Group: sid (26 Levels)
    Name1                Name2                Type          Estimate    Lower       Upper   
    '(Intercept)'        '(Intercept)'        'std'           10.286      7.3305      14.433
    'T'                  '(Intercept)'        'corr'         -0.9923    -0.99685    -0.98125
    'siteN'              '(Intercept)'        'corr'         0.42167    -0.13976     0.77794
    'T'                  'T'                  'std'          0.22742     0.16139     0.32047
    'siteN'              'T'                  'corr'        -0.44807    -0.79444     0.11823
    'siteN'              'siteN'              'std'          0.10547    0.065159     0.17072

Group: runN (6 Levels)
    Name1                Name2                Type          Estimate    Lower        Upper   
    '(Intercept)'        '(Intercept)'        'std'           1.8846      0.94663      3.7521
    'T'                  '(Intercept)'        'corr'        -0.99663          NaN         NaN
    'siteN'              '(Intercept)'        'corr'         0.98096      0.98027     0.98163
    'T'                  'T'                  'std'          0.04452      0.02311    0.085765
    'siteN'              'T'                  'corr'        -0.96172          NaN         NaN
    'siteN'              'siteN'              'std'         0.030653    0.0068465     0.13724

Group: Error
    Name             Estimate    Lower     Upper 
    'Res Std'        1.3459      1.2898    1.4044

%}