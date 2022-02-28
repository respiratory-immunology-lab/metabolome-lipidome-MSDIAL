# Secondary MS1 Lipid identification with Lipidmaps
After processing the Lipidomic LCMS data using the MS-DIAL pipeline, we can obtain secondary annotations for MS1 data by matching the mass data to LIPID MAPS Structure Database (LMSD). Before this step, you should run your data through the GNPS MS/MS secondary annotation steps (workflow here).

The LIPID MAPS Structure Database (file version LMSD 2022-02-16) was formated and saved as a .RDS file ([available here](https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/blob/main/LIPIDMAPS_processing/LMSD_2022_02_16.Rds)).
If a new version of the LMSD is available you can find the code for formatting [here](...).

To facilitate adding these annotations to the SummarizedExperiment object, a function is available at the top of the page. Its use will be explained here.

# Method

You should currently have a SummarizedExperiment object that has been pre-processed with pmp. The next stage is to match the MS1 mass data for each feature to the LMSD database file.

Appending LMSD annotations to SummarizedExperiment objects.

The first step is to load the formatted LMSD data.frame into your R session.

```{r}
# Load HMDB dataset
lmsd_df <- readRDS(here::here('lmsd', 'LMSD_2022_02_16.Rds.rds'))

### OR

lmsd_df <- readRDS('/path/to/lmsd/dataframe/')
```

Then, using the add_lipidmaps() function available [here]() we can search LMSD annotations in the data.frame and add them to our SummarizedExperiment objects, as shown below with an example BAL lipidomics object.


If you are using the output of the pmp_preprocess() function, you should annotate both the glog_results and imputed_results list components for consistency in case you need to use the imputed values downstream instead of the glog-transformed data.

The function takes three parameters:

* `metab_SE`: the `SummarizedExperiment` object.
* `lmsd`: the formatted LMSD data.frame.
* `mass_tol`: the MS1 mass tolerance value (set to 0.002 Da by default).
* `cores`: uses `doParallel` and `foreach`, either select number of cores or leave as NA and will use all cores-1. 

```
# Search annotations in LMSD and add to the SE objects
# Create list for all distinct Lipid Maps matching mz in tolerance range 0.002, an aggregated df of distinct lipids and a df to replace SummarizedExperiment metadata [rowData(metab_glog)]
lmsd_ann_list <- add_lmsd(metab_SE=metab_glog, lmsd=lmsd_df, mass_tol = 0.002) 

# use metadata_lmsd_table to replace SE object metadata
rowData(metab_glog) <- lmsd_ann_list$metadata_lmsd_table
```

# Comparing annotations from MS-DIAL, LMSD, GNPS, HMDB, and KEGG

To compare the assigned annotations from each of the methods the function `compare_annotations()`, available in this [folder](). It will produce a data.frame containing only features with at least one annotation, and allow us see whether the annotations typically agree with each other.

The function takes two argument: 
* `metab_SE`: a SummarizedExperiment object that has undergone secondary annotation with LMSD, GNPS and HMDB. 
* `agg_lmsd_ann_table`: the aggregated data.frame of distinct lipids create by `add_lmsd()`

```
# Prepare data.frame with alignment IDs and all annotations and filter for at least one annotation
msdial_lmsd_gnps_hmdb <- compare_annotations(metab_SE = metab_glog, 
                                             agg_lmsd_ann = lmsd_ann_list$agg_lmsd_df)
```

An example of the output should looks like:

| Command | Description |
| --- | --- |
| git status | List all new or modified files |
| git diff | Show file differences that haven't been staged |



