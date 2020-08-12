% This script should be updated and rerun whenever a new dataset is added to the repo.
% append a name to st_dataset variable, add a section to model your dataset
% dataset below following the template provided as a comment block, 
%
% This script critically depends on the <datasetname>_dataset_obj.mat files
% in datasets/behavior, not on the fMRI data, so make sure you have your
% *dataset_obj.mat file in the right place. You should not be adding
% datasets to this repo without a dataset_obj.mat file.

%
% We then compute SS_regressor following the method for computing R2
% detailed here:
% https://www.researchgate.net/publication/306347340_A_Pseudo_Decomposition_of_R2_in_Multiple_Linear_Regression
% and here:
% http://biol09.biol.umontreal.ca/borcardd/partialr2.pdf (no. 1)

close all;
clear all;

warning('off','all')

addpath(genpath('/projects/bope9760/CanlabCore/CanlabCore'));
addpath(genpath('/projects/bope9760/spm12'));
addpath(genpath('/projects/bope9760/single_trials_overview/repo')); % single_trials repo
addpath('/work/ics/data/projects/wagerlab/labdata/projects/canlab_single_trials_for_git_repo/'); % single trial data on blanca

%% import and prep data
% - concatenate subject data matrices from canlab_dataset objects
%
% - remove trials with high_vif or vif > 2.5, for consistency with subsequent
%   imaging datasets
%
% - separate pain from non-pain datasets
%
% - remove trials with missing responses (rating == nan)
%
% - add subject and study ids

st_datasets = {'nsf','bmrk3','bmrk4','scebl','ie2','ie',...
    'exp','stephan','romantic','ilcp'};
n_datasets = length(st_datasets);

dat = cell(n_datasets,1);
for i = 1:length(st_datasets)
    this_dat = importdata([st_datasets{i}, '_dataset_obj.mat']);
    dat{i} = array2table(cell2mat(this_dat.Event_Level.data'),...
        'VariableNames',this_dat.Event_Level.names);
    uniq_subject_id = this_dat.Subj_Level.id;
    subject_id = [];
    for j = 1:length(uniq_subject_id)
        subject_id = [subject_id, repmat(uniq_subject_id(j),1,size(this_dat.Event_Level.data{j},1))];
    end
    dat{i}.subject_id = subject_id';
    % this removes non-thermal data, we'll add these back in later
    dat{i} = dat{i}(~isnan(dat{i}.rating) & ~isnan(dat{i}.T),:);
    
    if contains(st_datasets{i},'bmrk5')
        dat{i}.soundintensity = [];
    end
    fnames = dat{i}.Properties.VariableNames;
    if any(ismember(fnames,'high_vif'))
        dat{i}(dat{i}.high_vif == 1,:) = [];
    end
    if any(ismember(fnames,'vif'))
        dat{i}(dat{i}.vif >= 2.5,:) = [];
    end
end

% remove warm (non-painful) trials (rating <= 100)
% note, this makes it harder to interpret some self regulation effects,
% since many regulate down conditions end up getting dropped this way, but
% because outcome measures mean different things at sub 100 vs. over 100,
% we drop them for linear outcome modeling purposes
bmrk3rating = dat{ismember(st_datasets,'bmrk3')}.rating - 100;
dat{ismember(st_datasets,'bmrk3')}.rating = bmrk3rating;
dat{ismember(st_datasets,'bmrk3')} = dat{ismember(st_datasets,'bmrk3')}(bmrk3rating > 0,:);

% take care of non-pain data
this_dat = importdata('bmrk3_dataset_obj');
dat{end+1} = array2table(cell2mat(this_dat.Event_Level.data'),...
        'VariableNames',this_dat.Event_Level.names);
uniq_subject_id = this_dat.Subj_Level.id;
subject_id = [];
for j = 1:length(uniq_subject_id)
    subject_id = [subject_id, repmat(uniq_subject_id(j),1,size(this_dat.Event_Level.data{j},1))];
end
dat{end}.subject_id = subject_id';
dat{end} = dat{end}(dat{end}.rating <= 100,:);

this_dat = importdata('bmrk5_dataset_obj');
dat{end+1} = array2table(cell2mat(this_dat.Event_Level.data'),...
        'VariableNames',this_dat.Event_Level.names);
uniq_subject_id = this_dat.Subj_Level.id;
subject_id = [];
for j = 1:length(uniq_subject_id)
    subject_id = [subject_id, repmat(uniq_subject_id(j),1,size(this_dat.Event_Level.data{j},1))];
end
dat{end}.subject_id = subject_id';
dat{end} = dat{end}(isnan(dat{end}.T) & ~isnan(dat{end}.rating),:);
dat{end}.T = [];

%c = seaborn_colors;
%c = shuffles(c);

%% plot experimental variance in each dataset
% - We model subjects as fixed effects
%
% - We use backwards stepwise regression to select additional experimental 
%   factors
%
% - Not all subject have trial sequence data, but for datasets where this
%   is missing we provide the mean
%
% - model interactions in cases where these were key factors (e.g.
%   plaebo*value in stephan's placebo dataset), but otherwise stick to
%   linear effects.
%
% - adjust labels in pie plots to avoid overlapping numbers
%
% - use consistent units across datasets, and concatenate all data for a
%   final aggregate analyses (top left pie plot). Use within-dataset
%   z-scored pain ratings for aggregate analyses. Exclude non pain data.
%
% - standardize all variables before regression modeling to get
%   standardized regression coefficients (necessary for partial R2
%   calculation)
%
% Note that for agregate analysis "subject fixed effects" are colinear with
% study effects, so you should expect the "subject effects" to explain a
% disproportionate amount of the variance relative to what they explain in 
% individual dataset analyses.
c = colormap('lines');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TEMPLATE (copy for new datasets) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
this_st_id = ???; % index corresponding to entry in st_dataset variable
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);
IV = ???; % cell array containing column header names of t to treat as IV

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar', IV,...
    'Intercept',false);

% compute partial r2
partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

% here we create a pie graph with four slices, experimental stimulus 
% components, sensitization and habituation components, cognitive
% manipulation components and subject effects. These are just different
% sums of R2 from the R2 variable above. This is mostly done automatically
% below, but check the variable names in m to ensure all relevant variables
% are accounted for below. you may need to modify this slightly.
figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:4
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = 'white'; % Subj fixed fx variance
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));
%}

%
% eval nsf
%
this_st_id = 1;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','overallN','siteN','runN'})),...
    'Intercept',false);

% compute partial r2
partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:4
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = 'white'; % Subj fixed fx variance
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% eval bmrk3pain
%
this_st_id = 2;
newdat_st = dat{this_st_id};

newdat_st = newdat_st(~newdat_st.high_vif,:);
[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','siteN','overallN','high_vif'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding','reg'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(3,1));
for i = 2:2:6
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(3).FaceColor = 'white'; % Subj fixed fx variance
f(2).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), 'pain\n'], '\n']));

%
% eval bmrk4
%
this_st_id = 3;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','siteN','runN','high_vif'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:8
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f(4).Position = f(4).Position + [0,0.25,0];
f(6).Position = f(6).Position - [0,0.25,0];
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(4).FaceColor = 'white'; % Subj fixed fx variance
f(3).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% bmrk5pain
%
this_st_id = 4;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','vif'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:6
    if strmatch(f(i).String,'0%')
        f(i).String = '';
    end
end
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(3).FaceColor = 'white'; % Subj fixed fx variance
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), 'pain\n'], '\n']));

%
% remi
%
this_st_id = 5;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','runN','overallN','vif','placebo'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
drug_idx = ismember(partialR2.Row,{'drug'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding','open'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx)),sum(R2(drug_idx))],ones(5,1));
for i = 2:2:10
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f(4).Position = f(4).Position - [0.2,0,0];
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(4).FaceColor = 'white'; % Subj fixed fx variance
f(3).FaceColor = c(1,:); % cog variance 
f(5).FaceColor = c(2,:); % drug variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% eval scebl
%
this_st_id = 6;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','overallN','cue','vif'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:8
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f(4).Position = f(4).Position - [0.25,-0.25,0];
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(4).FaceColor = 'white'; % Subj fixed fx variance
f(3).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% eval ie2
%
this_st_id = 7;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','high_vif','overallN'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:8
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f(4).Position = f(4).Position - [0.25,-0.25,0];
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(4).FaceColor = 'white'; % Subj fixed fx variance
f(3).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% eval ie
%
this_st_id = 8;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','siteN','runN','overallN'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:6
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(3).FaceColor = 'white'; % Subj fixed fx variance
f(2).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% eval_exp
%
this_st_id = 9;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','runN','siteN','vif'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:8
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(4).FaceColor = 'white'; % Subj fixed fx variance
f(3).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% eval levoderm
%
this_st_id = 10;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','vif','reveal','conditioningN','runN','siteN','overallN'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding','reveal'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:6
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(3).FaceColor = 'white'; % Subj fixed fx variance
f(2).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% eval stephan
%
this_st_id = 11;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

t.int = zscore(t.placebo.*t.value);
cov_names = [cov_names,'int'];

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'T','sid','rating','overallN','vif','runN'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','int','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:6
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f = f(1:2:end);
f(1).FaceColor = c(3,:); % Sensitization/Habituation
f(3).FaceColor = 'white'; % Subj fixed fx variance
f(2).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% eval romantic
%
this_st_id = 12;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','overallN','runN','T','high_vif'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:4
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f(2).Position = f(2).Position - [-0.25,0,0];
f(4).Position = f(4).Position - [0.25,0,0];
f = f(1:2:end);
f(1).FaceColor = c(3,:); % Sensitization/Habituation
f(3).FaceColor = 'white'; % Subj fixed fx variance
f(2).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% eval ilcp
%
this_st_id = 13;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','cue','overallN','runN','vif'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:8
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(4).FaceColor = 'white'; % Subj fixed fx variance
f(3).FaceColor = c(1,:); % cog variance 
title(sprintf([[strrep(strrep(st_datasets{this_st_id},'_data.mat',''),'_',' '), '\n'], '\n']));

%
% bmr3heat
%
this_st_id = 14;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','high_vif','overallN','siteN'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:6
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f(2).Position = f(2).Position - [-0.25,0,0];
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(3).FaceColor = 'white'; % Subj fixed fx variance
title(sprintf('bmrk3warm\n\n'));

%
% bmrk5snd
%
this_st_id = 15;
newdat_st = dat{this_st_id}; 

[t, cov_names] = get_std_tbl_with_fixed_fx(newdat_st);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,{'sid','rating','high_vif','vif','overallN','siteN'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T','soundintensity'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), this_st_id+1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx))],ones(4,1));
for i = 2:2:6
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f(2).Position = f(2).Position - [-0.3,0,0];
f(4).Position = f(4).Position - [0.3,0,0];
f = f(1:2:end);
f(1).FaceColor = c(4,:); % sound
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(3).FaceColor = 'white'; % Subj fixed fx variance
title(sprintf('bmrk5snd\n\n'));

%
% eval all pain
%
t = expand_metadata_table(dat{1:length(st_datasets)});
for i = 1:length(t)
    t{i}.st_id = cellstr(repmat(st_datasets{i},height(t{i}),1));
    t{i}.rating = zscore(t{i}.rating);
end
t = tbl_vcat(t{:});
[~,~,sid] = unique([char(t.subject_id), char(t.st_id)],'rows','stable');
t.subject_id = sid;
[~,~,st_id] = unique(char(t.st_id),'rows','stable');
t.st_id = st_id;

% manually set some defaults
t.placebo(isnan(t.placebo)) = 0;
t.drug(isnan(t.drug)) = 0;
t.open(isnan(t.open)) = 0;
t.open(t.placebo == 1) = 1; % all placebos are effectively "open" drug label in the sense used in the remifentanil study
t.handholding(isnan(t.handholding)) = 0;
t.reg(isnan(t.reg)) = 0;
t.reveal(isnan(t.reveal)) = 0;
t.social(isnan(t.social)) = 0;
t.value(isnan(t.value)) = 0;
t.ctrl(isnan(t.ctrl)) = 0;

t.high_vif(isnan(t.high_vif)) = 0;
t.vif(isnan(t.vif)) = 0;

new_t = t(~t.high_vif | t.vif > 2.5,:);
t = new_t;

t.int = t.value.*t.placebo;

[t, cov_names] = get_std_tbl_with_fixed_fx(t);

m = fitlm(t,'ResponseVar','rating',...
    'PredictorVar',cov_names(~ismember(cov_names,...
    {'conditioningN','high_vif','rating',...
    'soundintensity','st_id','subject_id','vif','reveal'})),...
    'Intercept',false);

partialR2 = partialRsquared(m);
if abs(sum(partialR2.partialR2) - m.Rsquared.Ordinary) > 0.0001
    warning(['Partial R2 computation doesn''t seem to be working for study ' int2str(this_st_id)]);
end

R2 = table2array(partialR2);

figure(1)
T_idx = ismember(partialR2.Row,{'T'});
sens_idx = ismember(partialR2.Row,{'siteN','overallN','runN'});
sid_idx = find(contains(partialR2.Row,'sub'));
cog_idx = ismember(partialR2.Row,{'cue','ctrl','social','placebo','value','reveal','handholding'});
drug_idx = ismember(partialR2.Row,'drug');
subplot(ceil(sqrt(1+n_datasets)), ceil((1+n_datasets)/ceil(sqrt(1+n_datasets))), 1);
f = pie([sum(R2(T_idx)),sum(R2(sens_idx)),sum(R2(cog_idx)),sum(R2(sid_idx)),sum(R2(drug_idx))],ones(5,1));
for i = 2:2:10
    if strmatch(f(i).String,'0%')
        f(i).String = ''
    end
end
f(4).Position = f(4).Position - [0.1,-0.1,0];
f(6).Position = f(6).Position - [0.1,0.2,0];
f = f(1:2:end);
f(1).FaceColor = c(7,:); % Temp
f(2).FaceColor = c(3,:); % Sensitization/Habituation
f(4).FaceColor = 'white'; % Subj fixed fx variance
f(3).FaceColor = c(1,:); % cog variance
f(5).FaceColor = c(2,:); % drug variance
title(sprintf('Pain Datasets (n=13) \n Aggregate\n'))

l = legend({'T','Sens/Habit','Cog','Subject','Drug'},'location','southwest');
l.Position = l.Position - [0.17,0,0,0];
set(get(gca,'Legend'),'FontSize',10)

p = get(gcf,'Position');
set(gcf,'Position',[p(1:2), 1248, 966]);

%% Quantify signal strength using SVR brain decoding
% - obfuscate dataset identity to discourage chery picking for future
%   analysis based on SNR. Any studies with SNR so low as to suggest
%   corrupt data will be identified. Otherwise, the purpose of what follows
%   is to provide guidelines for what one can expect in terms of decoding
%   performance when using linear MVPA methods.
%
% - control for subject fixed effects, since experience suggests between
%   subject effects aren't detected very well. We will print out raw r^2
%   values for predicted vs. observed though for comparison for those who
%   are curious.

st_datasets = {'nsf','bmrk3pain','bmrk4','bmrk5pain','remi','scebl','ie2','ie',...
    'exp','levoderm','stephan','romantic','ilcp','bmrk3warm','bmrk5snd'};

% obfuscate identity
st_datsets = st_datasets(randperm(length(st_datasets)));

[cverr, stats, optout, obs, subject_id] = deal(cell(length(st_datasets),1));
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'))
end
parpool(16);
parfor i = 1:length(st_datasets)
    warning('off','all')
    % capture output to prevent clues that might reveal which dataset this
    % is
    this_data = load_image_set(st_datasets{i});
    
    % remove any trials lacking responses (dat.Y = nan), or with high VIF
    good_trials = ones(length(this_data.Y),1);
    
    fnames = this_data.metadata_table.Properties.VariableNames;
    if any(ismember(fnames,'high_vif'))
        good_trials(this_data.metadata_table.high_vif == 1) = 0;
    end
    if any(ismember(fnames,'vif'))
        good_trials(this_data.metadata_table.vif >= 2.5) = 0;
    end
    good_trials(isnan(this_data.Y)) = 0;
    good_trials = logical(good_trials);
    
    this_data = this_data.get_wh_image(good_trials);
    
    % specify folds manually to maintain subject groupings across fold
    % slicings
    [~,~,subject_id{i}] = unique(char(this_data.metadata_table.subject_id),'rows','stable');
    cv = cvpartition2(ones(size(this_data.dat,2),1), 'KFOLD', 5, 'Stratify', subject_id{i});
    fold_labels = zeros(size(this_data.dat,2),1);
    for j = 1:cv.NumTestSets
        fold_labels(cv.test(j)) = j;
    end
    
    [cverr{i},stats{i},optout{i}] = this_data.predict('algorithm_name', 'cv_svr', 'nfolds', fold_labels, 'error_type', 'mse', 'useparallel', 0, 'verbose', 0);
    obs{i} = this_data.Y;
end

for i = 1:length(st_datasets)
    % plot results
    figure(2)
    subplot(ceil(sqrt(n_datasets)), ceil((n_datasets)/ceil(sqrt(n_datasets))), i);
    r2 = corr(stats{i}.yfit, obs{i}).^2;
    f = pie(r2);
    f = f(1:2:end);
    f(1).FaceColor = [0.5,0.5,0.5]; % svr
    if ismember(st_datasets{i},{'bmrk3warm','bmrk5pain'})
        title(sprintf('Random dataset %d\n(nonpain)\n', i));
    else
        title(sprintf('Random dataset %d\n\n', i));
    end
end
save('st_datasets.mat','st_datasets','cverr','stats','optout','obs','subject_id');

p = get(gcf,'Position');
set(gcf,'Position',[p(1:2), 1248, 966]);

%% Bad datasets
% -
% 
% -
