## Guidelines

1. Raw event level (single trial) beta maps only. No averaging over multiple events. 
2. Avoid adding new functions not directly related to data organization and management.
3. Keep the fmri_data_st class backwards compatible with fmri_data. No new properties.
4. No unpublished data. Include references in all fmri_data and canlab_data objects.
5. Use consistent metadata naming conventions. Check existing canlab_dataset objects to see if your metadata type is already included elsewhere, and if so use the same variable names.
6. fmri_data.dat should be linearly related to fmri_data.Y (to a first approximation). Either continuous outcomes from a single task or category labels for multiple types of tasks might satisfy this, but continuous outcomes of multiple tasks measures on the same scale (e.g. a pain intensity rating and sound intensity rating task, like in BMRK5) are unlikely to meet this requirement. Subdivide such data into two datasets, so that whatever you import with load_image_set() can have some sort of immediate usability without further preprocessing.
7. Document the creation of new datasets in the src folder so that if something crops up that's wrong, the source can be tracked down and fixed.
8. Include all data. If entries are missing anywhere (e.g. outcome measures or labes) use nans to indicate this.

## Instructions for adding data

Create a canlab_dataset object for your project. This should contain, if nothing else, event level (a.k.a. trial level) experimental manipulations, and a citation to a paper describing the experment. Put the data in datasets/behavior. src/datasets/behavior/TEMPLATE.m in this directory is provided to facilitate development of canlab_dataset objects from scratch. Keep your canlab_dataset creation script and upload it to this src directory too as a record of what was done. If updating existing canlab_dataset objects their original creation scripts are here to facilitate future additions.
- Include trial sequence information in your metadata for estimating sensitization and habituation effects. If unavailable fille overallN, siteN and runN with their mean values. siteN indicates trial indices. runN indicates block indices (e.g. if multiple sites are stimulated).

Upload your fri_data object to a public URL accessible repository. 
- Encrypt the file if privacy is desired using the AES class in opt/matlab_aes. For simplicity you can just use the lib/encrypt_dataset which uses a default repository wide password for encryption. Otherwise invoke the AES class with a password of your choosing. See 'help AES' for usage examples.
- Include all event level metadata in the metadata_table field (will be redundant with canlab_datasets).

create a load_*datasetName*.m file in datasets/brain/. This will get called now if you run canlabCore/load_image_set('*datasetName*'). Instruct it to searhch your matlab path for your dataset file, and once loaded to cast it to fmri_data_st() before returning.
 - Add instructions for downloading the dataset if it's found missing. Alternatively invoke lib/download_dataset(*datasetName*) and add a case for your dataset therein. Require passing 'forcedl' to load_*datasetName*(), or a user confirmation before downloading to avoid surprising users with unexpected large files occupying limited disk space.
- If your uploaded dataset uses encryption, add instructions for decrypting your dataset. You can invoke lib/decrypt_dataset if you used lib/encrypt_dataset for encryption (since both assume the default password). Otherwise invoke the AES class directly using your own password/hash.
- Consider coding an entry that performs an MD5 check on your (decrypted) data in load_*datasetName*.m. You may find datasets/brain/wrongmd5() helpful for this. Simply add your md5 sum to datasets/brain/checksums.md5 for compatibility.
- update load_all_single_trials.m to include a reference to your load_*datasetName*.m script.

update QC/QC_main.m to include an entry for behavioral modelling of new datasets and for SVR modeling of its brain data.
