% This script assembles fmri_data objects containing the brain info from
% the single trial dataset google drive (mostly, bmrk5 is the exception,
% but it's trivial one, as the same data is on the google drive), and takes
% what are multiple files there and combines them into a single dataset
% object.
%
% note the second argument of prep_dataset() is not used if the third
% argument is provided. The second argument is a study specific metadata
% file from Tor's google drive of single trial data. The third argument is
% the canlab_dataset object provided in this repo, so you should have that.
% if you can get the raw image files either from the google drive or
% elsewhere, you should be able to rerun these. Many (but not all) the
% files provided in the second argument are also available in the google
% drive.
%
% Update savePath
% Update study specific paths (e.g. nsfPath, bmrk3Path, etc) to point to
% where your single trial image files are
%
% update library import statements to point to the dependencies locations
% on your system

clear all; close all;

%% update these
analysisRoot    = '/projects/bope9760/single_trials_overview/';
savePath        = [analysisRoot,'resources/'];
addpath(genpath('/projects/bope9760/CanlabCore/CanlabCore'));
addpath(genpath('/projects/bope9760/spm12'));
addpath(genpath([analysisRoot,'/repo'])); % single_trials repo

nsfPath             = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/nsf/';
bmrk3Path           = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/Tor_bmrk3_datashare/';
bmrk4SmoothPath     = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/bmrk4_smoothed_withbasis/';
remiPath            = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/REMI_wani/';
ilcpPath            = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/ILCP_wani/';
sceblPath           = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/SCEBL_single_trial_Leonie/';
ie2Path             = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/IE2_NEW/';
iePath              = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/ie_for_tor/';
expPath             = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/Expectancy/';
levodermPath        = '/work/ics/data/projects/wagerlab/labdata/current/bogdan_spillover_storage/Wagerlab_Single_Trial_Pain_Datasets/Data/levoderm/';
stephanPath         = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/placebo_value_stephan/';
romanticPainPath    = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/romantic_pain/';
bmrk5Path           = '/work/ics/data/projects/wagerlab/labdata/projects/bmrk5_painsound/';

% fixes CanlabCore & spm12 compatibility issue
try
    rmpath(genpath(spider_path));
end

%% prep NSF data (NEJM and Atlas 2014 Pain paper)
% pain scale may only be 0-8 in this one. See p2 of Atlas 2015 Pain.
fprintf('Preparing NSF data...\n');

tmpData = prep_dataset(nsfPath, 'meta.mat',importdata('nsf_dataset_obj.mat'));
tmpData.source_notes = 'NSF data aggregated from Tor Wager''s single trials Google Drive';
tmpData.Y_descrip = 'Pain';

for i = 1:size(tmpData.dat,2)
    tmpData.dat(isnan(tmpData.dat(:,i)),i) = 0;
end
tmpData = tmpData.remove_empty;

save([savePath,'nsf_data.mat'],'tmpData','-v7.3');
clear tmpData;

%% prep BMRK3 data (NEJM paper, but also Wani's 2015 PLoS paper)
fprintf('Preparing BMRK3 data...\n');
% ratings below 100 are nonpainful

tmpData = prep_dataset(bmrk3Path, 'bmrk3_single_trial_model.mat',importdata('bmrk3_dataset_obj.mat'));
if any(isnan(tmpData.Y))
    warning(['Dropping ' sum(isnan(tmpData.Y)) ' trials from data due to failure to disambiguate heat from pain']);
end
bmrk3pain = fmri_data(fmri_data_st(tmpData).get_wh_image(find(tmpData.Y > 100)));
bmrk3warm = fmri_data(fmri_data_st(tmpData).get_wh_image(find(tmpData.Y <= 100)));

bmrk3pain.source_notes = 'bmrk3pain img data from Tor Wager''s single trials Google Drive. Metadata also from wagerlab/labdata/current/BMRK3/ HPC storage';
bmrk3warm.source_notes = 'bmrk3warm img data from Tor Wager''s single trials Google Drive. Metadata also from wagerlab/labdata/current/BMRK3 HPC storage';

bmrk3pain.Y_descrip = 'Pain';
bmrk3warm.Y_descrip = 'Warmth';

if size(bmrk3pain.dat,2) + size(bmrk3warm.dat,2) ~= size(tmpData.dat,2)
    warning('Trials dropped from bmrk3');
end

%fmri_data/cat doesn't handle images_per_session field well.
% lvls = [training_data.images_per_session(:); tmpData.images_per_session(:)];
% training_data = training_data.cat(tmpData);
% training_data.images_per_session = lvls;
save([savePath,'bmrk3pain_data.mat'],'bmrk3pain','-v7.3');
save([savePath,'bmrk3warm_data.mat'],'bmrk3warm','-v7.3');

clear tmpData bmrk3pain bmrk3warm;

%% prep BMRK4 smoothed data (Krishnan 2016 eLife paper, somatic vs vicarious pain)
fprintf('Preparing BMRK4 data...\n');

tmpData = prep_dataset(bmrk4SmoothPath, 'bmrk4_single_trial_model.mat',importdata('bmrk4_dataset_obj.mat'));
tmpData.source_notes = 'BMRK4 img data from Tor Wager''s single trials Google Drive (bmrk4_smoothed_with_basis). Metadata also from wagerlab/labdata/current/bmrk4/Data/BMRK_Data HPC resources.';
tmpData.Y_descrip = 'Pain';

save([savePath,'bmrk4_data.mat'],'tmpData','-v7.3');

clear tmpData;

%% prep REMI data (Atlas 2012 J Neuro, remifentanil vs. expectancy paper)
% 0 - 8 VAS scale
fprintf('Preparing REMI data...\n');

tmpData = prep_dataset(remiPath, 'remi_st_metadata.mat',importdata('remi_dataset_obj.mat'));
tmpData.source_notes = 'Remi data from Tor Wager''s single trials Google Drive.';
tmpData.Y_descrip = 'Pain';

if isempty(tmpData.additional_info)
    tmpData.additional_info = struct('reference',{'Atlas, et al. (2012) Journal of Neuroscience'});
end


save([savePath,'remi_data.mat'],'tmpData','-v7.3');

clear tmpData;

%% prep APP-fMRI data
% I have no idea what this dataset is, and I don't have stimulus info for
% it, so skip...
%appPath        = '/projects/bope9760/Wagerlab_Single_Trial_Pain_Datasets/Data/APP-fMRI/';
%
%[tmpData, app_stim_cnt] = prep_dataset(appPath, 'mediation_variables2.mat');
%tmpData.Y = tmpData.Y;
% training_data = tmpData;
%save([savePath,'appfMRI_data.mat'],'tmpData','-v7.3');
%
%clear tmpData;
%% prep ILCP_wani data (Wani SIIPS paper, 2017 Nat Comm, study 6)
% (see Notes_on_studies.rtf in Wagerlab_Single_Trial_Pain_Datasets before
% considering other versions of the ILCP dataset.
fprintf('Preparing ILCP data...\n');

tmpData = prep_dataset(ilcpPath, 'ilcp_metainfo_waniupdate.mat',importdata('ilcp_dataset_obj.mat'));
tmpData.source_notes = 'ILCP data from ILCP_wani on Tor Wager single trials Google Drive.';
tmpData.Y_descrip = 'Pain';

save([savePath,'ilcp_data.mat'],'tmpData','-v7.3');

clear tmpData;
%% prep SCEBL data (Leoni unpublished)
% (see Notes_on_studies.rtf in Wagerlab_Single_Trial_Pain_Datasets before
% considering other versions of the ILCP dataset.
%
% Careful, these could be pain affect (0-5 scale), not pain intensity, at
% least iuf you go based on the paper
fprintf('Preparing SCEBL data...\n');

tmpData = prep_dataset(sceblPath, 'SCEBLdata_forTor_N26.mat',importdata('scebl_dataset_obj.mat'));
tmpData.source_notes = 'SCEBL data from Tor Wager''s single trials Google Drive.';
tmpData.Y_descrip = 'Pain';

save([savePath,'scebl_data.mat'],'tmpData','-v7.3');

clear tmpData;
%% prep IE2 pain data (???)
% (see Notes_on_studies.rtf in Wagerlab_Single_Trial_Pain_Datasets before
% considering other versions of the IE dataset.
%
% Be very cautious using this data. Everyone has 48/49C stim apparently,
% but everyone also has pain ratings that are too similar. This is
% suspisious.
fprintf('Preparing IE2 data...\n');

tmpData = prep_dataset(ie2Path, 'ie2_metadata.mat',importdata('ie2_dataset_obj.mat'));
tmpData.source_notes = 'ie2 img data from Tor Wager''s single trials Google Drive ("IE2_NEW"). Metadata from wagerlab/labdata/current/ie2/ HPC storage.';
tmpData.Y_descrip = 'Pain';

save([savePath,'ie2_data.mat'],'tmpData','-v7.3');

clear tmpData;

%% prep IE pain data (Roy 2014 Nat Neuro apparently)
fprintf('Preparing IE data...\n');

tmpData = prep_dataset(iePath, 'ie_model.mat',importdata('ie_dataset_obj.mat'));
tmpData.source_notes = 'ie img data from Tor Wager''s single trials Google Drive';

save([savePath,'ie_data.mat'],'tmpData','-v7.3');

clear tmpData;

%% prep EXP pain data  (Atlas 2010, J Neuro)
% 0-8 VAS scale

fprintf('Preparing EXP data...\n');

tmpData = prep_dataset(expPath, 'exp_meta.mat',importdata('exp_dataset_obj.mat'));
tmpData.source_notes = 'EXP data from Tor Wager''s single trials Google Drive. Metadata also from wagerlab/labdata/current/Expectancy/, especially the subfolder Behavioral/Sas_proc_mixed/, HPC storage.';
tmpData.Y_descrip = 'Pain';

save([savePath,'exp_data.mat'],'tmpData','-v7.3');

clear tmpData;

%% prep levoderm pain data  (possibly Schafer 2015)
fprintf('Preparing levoderm data...\n');

tmpData = prep_dataset(levodermPath, 'levoderm.mat',importdata('levoderm_dataset_obj.mat'));
tmpData.source_notes = 'levoderm data from Tor Wager''s single trials Google Drive';
tmpData.Y_descrip = 'Pain';

save([savePath,'levoderm_data.mat'],'tmpData','-v7.3');

clear tmpData;

%% prep stephan pain data  (Geuter 2013 Neuroimage)
% this data only has one stimulus level, but has tremendous variability in
% pain (due to placebo). Could be very useful for analyses that don't
% average over stimulus level
fprintf('Preparing Stephan placebo data...\n');

tmpData = prep_dataset(stephanPath, 'placebo_value_meta_tor_bog_mod.mat',importdata('stephan_dataset_obj.mat'));
tmpData.source_notes = 'Stephan''s Placebo data from from Tor Wager''s single trials Google Drive';
tmpData.Y_descrip = 'Pain';

save([savePath,'stephan_data.mat'],'tmpData','-v7.3');

clear tmpData;
%% prep romantic pain data  (? marina paper ?)
fprintf('Preparing romantic pain data...\n');

tmpData = prep_dataset(romanticPainPath, 'romantic_pain.mat',importdata('romantic_dataset_obj.mat'));
tmpData.source_notes = 'Romantic Pain data from Tor Wager''s single trials Google Drive';
tmpData.Y_descrip = 'Pain';

save([savePath,'romantic_data.mat'],'tmpData','-v7.3');
clear tmpData;

%% prep bmrk5
evmeta = importdata([bmrk5Path, 'single_trial_output_scnlab_july2016/single_trial_bmrk5_meta.mat']);

nsubj = size(evmeta.beta_name,1);

subj = cell(nsubj,1);
if isempty(gcp('nocreate'))
    parpool(16);
end
load('bmrk5_dataset_obj.mat','dataset_obj');
parfor i = 1:nsubj
    trials = find(cellfun(@(x1)(~isempty(x1)),evmeta.beta_name(i,:)));
    ntrials = length(trials);
    img = cell(ntrials,1);
    for j = 1:ntrials
        this_trial = trials(j);
        this_beta = evmeta.beta_name{i,this_trial};
        fname = [bmrk5Path, 'first_level_all_trs_single_trials/', this_beta];
        img{j} = fmri_data(fname);
    end
    subj{i} = cat(img{:});
    
    ad_info = dataset_obj.Event_Level.data{i};
    % include all ad_info that have entries (many will be
    % empty columns because we have columns for any factor from any
    % study, and no study has all factors)
    has_ad_info = any(cell2mat(dataset_obj.Event_Level.data'));
    t1 = array2table(ad_info(:,has_ad_info),...
        'VariableNames', dataset_obj.Event_Level.names(has_ad_info));
    t2 = table(cellstr(repmat(evmeta.subjects{i},ntrials,1)),'VariableNames',{'subject_id'});
    t = horzcat(t1,t2);
    
    subj{i}.Y = t.rating;    
    subj{i}.metadata_table = t;
    subj{i} = fmri_data_st(subj{i});
end

bmrk5 = cat(subj{:});
bmrk5.additional_info = struct('references',dataset_obj.Description.references);

bmrk5pain = fmri_data(bmrk5.get_wh_image(find(~isnan(bmrk5.metadata_table.T))));
bmrk5snd = fmri_data(bmrk5.get_wh_image(find(~isnan(bmrk5.metadata_table.soundintensity))));

bmrk5pain.source_notes = 'bmrk5pain data from wagerlab/labdata/projects/bmrk5_painsound/. Img data from the first_level_all_trs_single_trials/ subjfolder. Metadata from the single_trial_output_scnlab_july2016/ subfolder';
bmrk5snd.source_notes = 'bmrk5snd data from wagerlab/labdata/projects/bmrk5_painsound/. Img data from the first_level_all_trs_single_trials/ subjfolder. Metadata from the single_trial_output_scnlab_july2016/ subjfolder';

bmrk5pain.Y_descrip = 'Pain';
bmrk5snd.Y_descrip = 'Pain';

save([savePath,'bmrk5pain_data.mat'],'bmrk5pain','-v7.3');
save([savePath,'bmrk5snd_data.mat'],'bmrk5snd','-v7.3');
