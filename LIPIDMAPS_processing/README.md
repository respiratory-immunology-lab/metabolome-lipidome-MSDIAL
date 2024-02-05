# Secondary MS1 Lipid identification with Lipidmaps
After processing the Lipidomic LCMS data using the MS-DIAL pipeline, we can obtain secondary annotations for MS1 data by matching the mass data to LIPID MAPS Structure Database (LMSD). Before this step, you should run your data through the GNPS MS/MS secondary annotation steps (workflow here).

The LIPID MAPS Structure Database (file version LMSD 2022-02-16) was formated and saved as a .RDS file ([available here](./LMDB_231107.rds)).
If a new version of the LMSD is available you can find the code for formatting [here (still need to add)](...).

To facilitate adding these annotations to the SummarizedExperiment object, a function is available at the top of the page. Its use will be explained here.

# Method

You should currently have a SummarizedExperiment object that has been pre-processed with pmp. The next stage is to match the MS1 mass data for each feature to the LMSD database file.

Appending LMSD annotations to SummarizedExperiment objects.

The first step is to load the formatted LMSD data.frame into your R session.

```r
# Load LMSD dataset
lmsd_df <- readRDS(here::here('lmsd', 'LMSD_2022_02_16.Rds.rds'))

### OR

lmsd_df <- readRDS('/path/to/lmsd/dataframe/')
```

Then, using the add_lmsd() function available [here](./add_lmsd.R) we can search LMSD annotations in the data.frame and add them to our SummarizedExperiment objects, as shown below with an example BAL lipidomics object.


If you are using the output of the pmp_preprocess() function, you should annotate both the glog_results and imputed_results list components for consistency in case you need to use the imputed values downstream instead of the glog-transformed data.

The function takes three required parameters:

- `metab_SE`: the `SummarizedExperiment` object.
- `lmsd`: the formatted LMSD data.frame.
- `mass_tol`: the MS1 mass tolerance value (set to 0.002 Da by default).

To speed up the process, multiple cores can be used:

- `cores`: the number of parallel processes to run. Leave as NA and will use all cores-1.

```r
# Search annotations in LMSD and add to the SE objects
# Create list for all distinct Lipid Maps matching mz in tolerance range 0.002, an aggregated df of distinct lipids and a df to replace SummarizedExperiment metadata [rowData(metab_glog)]
lmsd_ann_list <- add_lmsd(metab_SE = metab_glog, 
                          lmsd = lmsd_df, 
                          mass_tol = 0.002,
                          cores = 8) 

# use metadata_lmsd_table to replace SE object metadata
rowData(metab_glog) <- lmsd_ann_list$metadata_lmsd_table
```

# Comparing annotations from MS-DIAL, LMSD, GNPS, HMDB, and KEGG

To compare the assigned annotations from each of the methods the function `compare_annotations_lip()`, available [here](./compare_annotations_lip.R). It will produce a data.frame containing only features with at least one annotation, and allow us see whether the annotations typically agree with each other.

The function takes two argument: 
* `metab_SE`: a SummarizedExperiment object that has undergone secondary annotation with LMSD, GNPS and HMDB. 
* `agg_lmsd_ann_table`: the aggregated data.frame of distinct lipids create by `add_lmsd()`

```r
# Prepare data.frame with alignment IDs and all annotations and filter for at least one annotation
msdial_lmsd_gnps_hmdb <- compare_annotations(metab_SE = metab_glog, 
                                             agg_lmsd_ann = lmsd_ann_list$agg_lmsd_df)
```
  
# Keeping only annotated features

From here, we can filter our SummarizedExperiment object for features with at least one annotation using the `keep_annotated_lip()` function [here](./keep_annotated_lip.R). It assigns rownames based on a naming hierarchy: LMSD > MS-DIAL > GNPS. A new rowData element will also be added to the `SummarizedExperiment` object called shortname.

This function will be applied to only a single SummarizedExperiment object, and the output can be assigned to a new object that can be used for downstream analyses. In the example below, we will run the function for just the glog-transformed data, and assign it to a new object called metab_glog.

The function takes only a single argument: a SummarizedExperiment object that has undergone secondary annotation with LMSD, GNPS and HMDB.

# Keep only annotated rows and generate shortname column
metab_glog <- keep_annotated(metab_glog$glog_results)
Aside from using rownames(), we can also retrieve the resulting "preferred" annotation names for our features as follows:

# Retrieve the shortname values
metab_shortnames <- rowData(metab_stool_glog)$shortname
While the shorter names are succinct and useful plotting, you can view the additional annotations at any time and alter as required.

# Get HMDB and KEGG annotations
hmdb_annotations <- rowData(metab_stool_glog)$HMDB
kegg_annotations <- rowData(metab_stool_glog)$KEGG

# Rights

Copyright (c) 2024 Respiratory Immunology lab, Monash University, Melbourne, Australia.
LMSD: LIPID MAPS® structure database. Sud M., Fahy E., Cotter D., Brown A., Dennis E., Glass C., Murphy R., Raetz C., Russell D., and Subramaniam S., Nucleic Acids Research 35, D527-32 (2006)
License: This pipeline is provided under the MIT license (See LICENSE.txt for details)
Authors: A. Butler, C. Pattaroni and M. Macowan






