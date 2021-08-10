# Secondary MS1 feature identification with HMDB

After processing our metabolomics LCMS data using the MS-DIAL pipeline, we can get secondary annotations for our MS1 data by matching the mass data to the [Human Metabolome Database (HMDB)](https://hmdb.ca/). Before this step, you should run your data through the GNPS MS/MS secondary annotation steps ([workflow here](https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/tree/main/gnps_processing)).

The HMDB database (version 4 - July 2021) has been formatted and filtered for features that are annotated or documented, and saved as a .RDS file ([available here](https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/blob/main/hmdb_processing/hmdb_detected_quantified_v4_20210701.rds)).

To facilitate the adding of these annotations to your `SummarizedExperiment` object, a function is available at the top of the page. Its usage will be explained here.

## Method

You should currently have a `SummarizedExperiment` object that has been pre-processed with [`pmp`](https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/tree/main/pmp_preprocessing), and then had the MS/MS data annotated with [GNPS](https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/tree/main/gnps_processing). The next stage is to match the MS1 mass data for each feature to the HMDB database file.

### Appending HMDB annotations to SummarizedExperiment objects

The first step is to load the formatted HMDB data.frame into your R session.

```R
# Load HMDB dataset
hmdb_df <- readRDS(here::here('hmdb', 'hmdb_detected_quantified_v4_20210701.rds'))

### OR

hmdb_df <- readRDS('/path/to/hmdb/dataframe/')
```

Then, using the `add_hmdb()` function available above, you can search the HMDB annotations in the data.frame and add them to your `SummarizedExperiment` objects, as shown below with an example stool metabolomics object.

If you are using the output of the `pmp_preprocess()` function, you should annotate both the `glog_results` and `imputed_results` list components for consistency in case you need to use the imputed values downstream instead of the glog-transformed data.

```R
# Search annotations in HMDB and add to the SE objects
metab_stool_pmp$glog_results <- add_hmdb(SEexperiment = metab_stool_pmp$glog_results,
                                         hmdb = hmdb_df, mass_tol = 0.002)
metab_stool_pmp$imputed_results <- add_hmdb(SEexperiment = metab_stool_pmp$imputed_results,
                                            hmdb = hmdb_df, mass_tol = 0.002)
```

### Comparing annotations from MS-DIAL, GNPS, HMDB, and KEGG



## Rights

* Copyright (c) 2021 Respiratory Immunology lab, Monash University, Melbourne, Australia.
* HMDB version 4.0: Wishart DS, Feunang YD, Marcu A, Guo AC, Liang K, et al., HMDB 4.0 — The Human Metabolome Database for 2018. Nucleic Acids Res. 2018. Jan 4;46(D1):D608-17. 29140435
* License: This pipeline is provided under the MIT license (See LICENSE.txt for details)
* Authors: C. Pattaroni and M. Macowan
