close all; clear all;

addpath(genpath('/projects/bope9760/single_trials_overview/repo')); % single_trials repo
addpath('/work/ics/data/projects/wagerlab/labdata/projects/canlab_single_trials_for_git_repo/'); % single trial data location

outPath = '/projects/bope9760/single_trials_overview/resources/'; % where to save encrypted files

st_datasets = {'bmrk3warm','bmrk5pain',...
    'bmrk5snd','remi','levoderm'};

for i = 1:length(st_datasets)
    this_ds = st_datasets{i};
    encrypt_dataset(which([this_ds,'_data.mat']), [outPath, this_ds, '.mat_encrypted']);
end
