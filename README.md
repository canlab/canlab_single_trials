# CANLab Single Trials Repository 

This repository is designed for education and model and algorithm development for multivariate pattern analysis (MVPA). Its purpose is to provide access to a currated dataset in a homogenous format suitable for immediate MVPA applications.

This repository will download (and decrypt) single trial data from the cloud if it is unavailable, and directly provides canlab_dataset objects with study metadata.

This repo also provides a set of convenience functions for loading and working with the single trials datasets (see "notable functions"), including several overloaded methods and classes designed as drop in replacements that will seemlessly integrate into the typical canlabCore workflow.

A simple demonstration of using this repo for MVPA algorithm development is illustrated in the [CANLab walkthrough script](https://canlab.github.io/_pages/canlab_single_trials_demo/demo_norming_comparison.html) at canlab.github.io.

## Data Overview
### Experimental factors affecting outcome measure (R<sup>2</sup>)
The figure below provides an overview of the effect of experimental factors on pain (or other outcome measures) for each dataset. Use this to select datasets appropriate to your applications (e.g. if you're interested in between subject effects, pick studies that show high subject variance). 

The models illustrated below estimate linear effects and subject means (fixed effects). Source [here](QC/QC_main.m). See the [published html](QC/html/QC_main.html) for details.

Caution: in the aggregated dataset analysis the "subjects" wedge includes dataset mean differences too since subjects and datasets are colinear. Do not naively assume it represents the mean of the individual dataset graphics.

![Image description](QC/html/QC_main_01.png)

## Setup 

Dependencies (add these to your path first, and ensure they're up to date)
- CanlabCore

Clone the repo, add it to your Matlab path. If single trial dataset are available, add those to your matlab path as well (otherwise they will be downloaded). 

If working on Discovery or Blanca please add the appropriate directory to your Matlab path rather than autodownloading, otherwise we'll end up with multiple redundant copies and waste shared resources:

- Discovery (Dartmouth): /dartfs/rc/lab/C/CANlab/labdata/projects/canlab_single_trials_for_git_repo
- Blanca (CU): /work/ics/data/projects/wagerlab/labdata/projects/canlab_single_trials_for_git_repo

Downloaded datasets come from the web in encrypted form and are automatically decrypted by downloader methods using built in credentials. Encryption capabilites are provided for developers, who should refer to the developer [guidelines](development.md).

## Usage

Example:

data = load_image_set('bmrk4', 'md5check')

See "Available Datasets" subsection for dataset names. 

### Notable functions ###

class fmri_data_st\
Inherets fmri_data object and overloads methods to better handle specific attributes of single trial datasets (e.g. metadata_table fields). You can use the datasets as fmri_data objects but you will have an easier time if they're cast as fmri_data_st instead (automatically done if called using load_image_set, see below))
- fmri_data_st/cat
- fmri_data_st/get_wh_image
- fmri_data_st/mean

function fmri_data_st/quantileByY()\
will return a dataset averaged over *.Y quantiles (e.g. quartiles). Useful for prototyping. Working with single trials directly is very slow and will sap your productivity. Don't do it.

class cvpartition2\
Inherits matlabs cvpartition, and uses the same invocation but respects block membership groupings (e.g. subject identity) when the stratify argument is specified. Useful for cross validated estimates on single trial data (where subjects should belong to either test or training splits, but not both).

extends load_image_set()
- Support importing single trial datasets by name. (e.g. nsf = load_image_set('nsf')). 
- Auto converts datasets to type fmri_data_st.
- Auto downloads missing datasets (after user prompt). 
- Optional: md5 check
- support importing all datasets in a single fmri_data_st object: dat = load_image_set('all_single_trials');

### Available fmri_data (or fmri_data_st) datasets 
The following are available as fmri_data objects (cast to fmri_data_st objects if imported using load_image_set()). Use explicitly with load_image_set(). e.g. load_image_set('nsf')
- nsf
- bmrk3pain
- bmrk3warm
- bmrk4
- bmrk5heat
- bmrk5snd
- exp
- ie
- ie2
- ilcp
- levoderm
- remi
- romantic
- scebl
- stephan

single trial fmri_data objects have,
- all trials (including non-response trials with nan entries), unmodified, taken from single trials google drive
- magnitude ratings (pain, heat or sound intensity for bmrk5) in fmri_data.Y field
- description of fmri_data.Y in fmri_data.Y_descrip
- trial metadata (e.g. cues, stimulus intensity, or other manipulations) in fmri_data.metadata_table field
- Relevant citations in fmri_data.additional_info.references
- source notes indicating where img and metadata were found (disambiguating where needed, e.g. bmrk4_smoothed_withbasis is indicated in source_notes for bmrk4).
- are subdivided by stimulus modality (e.g. bmrk3, bmrk5)

### Available canlab_dataset objects
\*Notice that dataset objects are not subdivided by modality.
- nsf
- bmrk3*
- bmrk4
- bmrk5*
- exp
- ie
- ie2
- ilcp
- levoderm
- remi
- romantic
- scebl
- stephan

dataset objects have,
- all relevant experimental variables, coded in a consistent fashion across studies
- coding scheme descriptions for any esoteric entries (e.g. 'open' in REMI or 'reveal' in stephan placebo)
- dataset references to cite
- source files for generating dataset objects, to trace any misentered data if any is found, or to facilitate addition of newfound data when it becomes available.
- for BMRK5, IE and IE2 there's age, gender, race and handedness information Subject Level information
- for REMI there's average dose Subject Level information
- all data was checked for sensible experimental effects using mixed models. Results documented as comments in src files

### Additional datasets
See [development](development.md) for guidelines regarding addition of new datasets

## Dataset sharing via github
While the purpose of this repository is not primarily data distribution, it does provide a model of how data can easily paired as a sidecar to a private github repository.

Datasets in this repository are not on github, but access to them is as-if they were. the load_\<datasetName\> functions automatically handle the missing data scenario by downloading encrypted datasets from a hardcoded location. Thus access to this repository is necessary and sufficient to acquire the data. 
  
This principle is easily extended to other datasets aside from those here, and I would encourage other people to adopt these methods and the code in this library to facilitate sharing companion data to github repositories, just don't put your data here if it's not consistent with the single trials format. 

We use a single encryption key for all datasets associated with this repository for simplicity, but the AES class provided here is compatible with arbitrary keys. Thus, if you want to adopt this repo's strategy for data sharing, simple take that class (and the javaMethodWrapper class dependency), encrypt a file with a password of your choice as indicated in the AES class help documentation, upload your encrypted file somewhere public, and create a matlab function (containing your password) which downloads and decrypts the file. Then any github repository with that file and the AES class becomes necessary and sufficient for accessing your data.

## Troubleshooting

If JavaMethodWrapper gives you an error when you try to load/download a dataset you should make sure you have the latest version of Matlab. Issues have been reported with Matlab 2015 under Windows and Matlab 2018b under mac, which were resolved with "in the AES.m I changed  import java.util.Base64; to  import java.util.Base64.\*; I think there was another one but I can't recall". This repo was developed under Matlab 2019b and redhat Linux.
