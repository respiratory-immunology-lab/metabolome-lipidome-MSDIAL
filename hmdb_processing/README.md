# Secondary MS1 feature identification with HMDB

After processing our metabolomics LCMS data using the MS-DIAL pipeline, we can get secondary annotations for our MS1 data by matching the mass data to the [Human Metabolome Database (HMDB)](https://hmdb.ca/). Before this step, you may decide to run your data through the GNPS MS/MS secondary annotation steps ([workflow here](https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/tree/main/gnps_processing)).

The HMDB database (version 5 - November 2023) has been formatted and filtered for features that are annotated or documented, and saved as an RDS file ([available here](./hmdb_metabolites_detect_quant_v5_20231102.rds)).

To facilitate the adding of these annotations to your `SummarizedExperiment` object, a function is available at the top of the page. Its usage will be explained here.

## Method

You should currently have a `SummarizedExperiment` object that has been pre-processed with [`pmp`](https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/tree/main/pmp_preprocessing), and then **optionally** had the MS/MS data annotated with [GNPS](https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/tree/main/gnps_processing). The next stage is to match the MS1 mass data for each feature to the HMDB database file.

### Appending HMDB annotations to SummarizedExperiment objects

The first step is to load the formatted HMDB data.frame into your R session.

```R
# Load HMDB dataset
hmdb_df <- readRDS(here::here('hmdb', 'hmdb_metabolites_detect_quant_v5_20231102.rds'))

### OR

hmdb_df <- readRDS('/path/to/hmdb/dataframe/')
```

Then, using the [`add_hmdb()` function](./add_hmdb.R), we can search the HMDB annotations in the data.frame and add them to our `SummarizedExperiment` objects, as shown below with an example stool metabolomics object.

If you are using the output of the `pmp_preprocess()` function, you should annotate both the `glog_results` and `imputed_results` list components for consistency in case you need to use the imputed values downstream instead of the glog-transformed data (although the standard practice is to use the glog-transformed data).

The function takes three required parameters:

- `metab_SE`: the `SummarizedExperiment` object.
- `hmdb`: the formatted HMDB data.frame.
- `mass_tol`: the MS1 mass tolerance value (set to 0.002 Da by default).

To speed up the process, multiple cores can be used:

- `cores`: the number of parallel processes to run.

```R
# Search annotations in HMDB and add to the SE objects
metab_stool_pmp$glog_results <- add_hmdb(metab_SE = metab_stool_pmp$glog_results,
                                         hmdb = hmdb_df, mass_tol = 0.002)
metab_stool_pmp$imputed_results <- add_hmdb(metab_SE = metab_stool_pmp$imputed_results,
                                            hmdb = hmdb_df, mass_tol = 0.002)
```

### Comparing annotations from MS-DIAL, GNPS, HMDB, and KEGG

We can now compare the assigned annotations from each of the methods using the function [`compare_annotations_met()`](./compare_annotations_met.R), available in this folder.
It will produce a data.frame containing only features with at least one annotation, and allow us see whether the annotations typically agree with each other.

The function takes only a single argument: a `SummarizedExperiment` object that has undergone secondary annotation with both GNPS and HMDB.

```R
# Prepare data.frame with alignment IDs and all four annotations and filter for at least one annotation
msdial_gnps_hmdb_glog <- compare_annotations_met(metab_stool_pmp$glog_results)
```

### Keeping only annotated features

From here, we can filter our `SummarizedExperiment` object for features with at least one annotation using the [`keep_annotated_met()` function](./keep_annotated_met.R).
It assigns rownames based on a naming hierarchy: HMDB > GNPS > MS-DIAL.
A new `rowData` element will also be added to the `SummarizedExperiment` object called `shortname`.

This function will be applied to only a single `SummarizedExperiment` object, and the output can be assigned to a new object that can be used for downstream analyses. In the example below, we will run the function for just the glog-transformed data, and assign it to a new object called `metab_stool_glog`.

The function takes only a single argument: a `SummarizedExperiment` object that has undergone secondary annotation with both GNPS and HMDB.

```R
# Keep only annotated rows and generate shortname column
metab_stool_glog <- keep_annotated(metab_stool_pmp$glog_results)
```

Aside from using `rownames()`, we can also retrieve the resulting "preferred" annotation names for our features as follows:

```R
# Retrieve the shortname values
metab_shortnames <- rowData(metab_stool_glog)$shortname
```

While the shorter names are succinct and useful plotting, you can view the additional annotations at any time and alter as required.

```R
# Get HMDB and KEGG annotations
hmdb_annotations <- rowData(metab_stool_glog)$HMDB
kegg_annotations <- rowData(metab_stool_glog)$KEGG
```

## Rights

* Copyright (c) 2024 Respiratory Immunology lab, Monash University, Melbourne, Australia.
* HMDB version 5.0: Wishart DS, Guo A, Oler E, Wang F, Anjum A, Peters H, Dizon R, Sayeeda Z, Tian S, Lee BL, Berjanskii M. HMDB 5.0: the human metabolome database for 2022. Nucleic acids research. 2022 Jan 7;50(D1):D622-31.
* License: This pipeline is provided under the MIT license (See LICENSE.txt for details)
* Authors: C. Pattaroni and M. Macowan
