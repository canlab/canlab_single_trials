clear all; close all;

analysisRoot    = '/projects/bope9760/single_trials_overview/';
addpath(analysisRoot)
init_script2
dataset_obj = canlab_dataset;

dataset_obj.Description.Experiment_Name = '';

fprintf(['Preparing ' dataset_obj.Description.Experiment_Name '  dataset obj...\n']);

%% update these paths as needed
dataRoot     = '';
% single trial image metadata
imgmeta        = importdata([dataRoot '/meta.mat']); 
% behavioral meta (might be identical to the above)
evmeta = importdata(''); 

%% ideally these will match across studies
dataset_obj.Subj_Level.names = {'age','male','race','right_handed'};
dataset_obj.Subj_Level.type = {'int','boolean','level','boolean'};
dataset_obj.Event_Level.names = {'rating', 'T', 'soundintensity', 'cue', ...
    'social', 'placebo', 'value', 'runN', 'siteN', 'overallN', 'vif', ...
    'ctrl', 'reveal','conditioningN','reg','handholding'};
dataset_obj.Event_Level.type = ...
    repmat('numeric',length(dataset_obj.Description.Event_Level));
dataset_obj.Event_Level.units = ...
    {'vas','C','levels','levels','levels','boolean','levels','count','count',...
    'count','unitless','boolean','boolean','count','levels','boolean'};

%% some study specific stuff
%datset_obj.Event_Level.descrip{ismember(dataset_obj.Event_Level.names,'high_vif')} = ...
%    'vif > 2.5';

%% get within subject data
nSubj = length(imgmeta.dat_obj);
hasData = zeros(nSubj,1);
subjVols = cell(nSubj,1);
volsPerSubj = zeros(nSubj,1);


for i=1:nSubj
    %% load subject's fmri_data obj
    subjNum = i;
    subj = imgmeta.dat_obj{subjNum};
    dataset_obj.Subj_Level.id = subj;
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
    these_cov = zeros( length(imgmeta.ratings{i}) ,length(dataset_obj.Event_Level.names));
    ...
    % assign covariate columns dynamically, not statically. If
    % dataset_obj.Event_Level.names changes you want this code here to
    % still work correctly. Changes to the column entries may occur as more
    % data is discovered or as changes are made to the dataset objects 
    % collection side. 
    %
    % additionally, code variables such that 0 corresponds to the default
    % in most studies. For instance, for placebo manipulations we code 0
    % for ctrl, and 1 for placebo, because in most studies there's no
    % placebo manipulation. This is helpful when merging studies with
    % different covariantes because it allows us to assign 0 as the default
    % value.
    
    dataset_obj.Event_Level.data{end+1} = these_cov;
end

if length(dataset_obj.Event_Level.data) ~= nSubj
    error('Length of dataset mismatch');
else
    save([analysisRoot 'resources/prep_canlab_dataset_objs/' ...
        dataset_obj.Description.Experiment_Name, '_dataset_obj.mat'],'dataset_obj');
end

%% sanity check
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

% center within subject since between subject factors are not as clear cut
% as within subject factors, so easier to interpret sensibility of results
% based on subject centered data (of course, exceptions should be made for
% a study design looking for between subject differences)
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
m = fitlme(tbl,'rating ~ T + siteN - 1 + (T + siteN - 1| sid) + (T + siteN | runN)','FitMethod','REML')

%{
    Save output from fitlme call here
%}