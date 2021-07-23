# Function to automate Peak Matrix Processing steps

The script in this folder will automate the initial pre-processing steps outlined in the `pmp_preprocessing` folder.

You can download the script `.R` file and place it in a folder within your R project.
To load the function into your R environment, you can then load it as follows:

```r
source(here::here('scripts', 'pmp_preprocess.R')

### OR

source('path/to/scripts/folder/pmp_preprocess.R')
```

### Dependencies

The following packages are required, and will be loaded by the function `pmp_preprocess()` if they are not currently loaded.

- `data.table`
- `tidyverse`
- `ggplot2`
- `pmp`
- `SummarizedExperiment`
- `S4Vectors`
- `ggsci`

### Usage

This function will take two output files from MS-DIAL and an associated sample metadata file (optional).
The two MS-DIAL output files required are the height tables from both the positive and negative ionisation modes.

The row names of the sample metadata file should match the column names of the samples in your MS-DIAL outputs only (not the blanks and QCs).
The function will reorder the metadata rows to match their order in the MS-DIAL outputs.

#### Function parameters

All of the parameters you can alter throughout the `pmp` pre-processing pipeline are available as arguments to the `pmp_preprocess()` function call.

Of particular note is the `samples_key` argument. 
If you do not manually indicate which columns contain your metabolite intensity values (using the `intens_cols` argument),
`pmp_preprocess()` will attempt to automatically determine the columns that contain your samples, blanks, and QCs.
To do this however, it requires that you give it something to search for in the column names: i.e. the prefix for your sample names.

Similarly, if you do not manually provide column indices containing the metabolite feature information (using the `info_cols` argument),
the function will automatically select the information columns using the 32 column names that MS-DIAL uses as of version 4.70.

The parameters you can set for the `pmp_preprocess()` function are as follows:

- `pos_df`: the MS-DIAL height output for the positive ionisation mode.
- `neg_df`: the MS-DIAL height output for the negative ionisation mode.
- `metadata`: you metadata `data.frame` for __samples__ only (e.g. patient information or clinical metadata for each sample)
- `samples_key`: a prefix used in the MS-DIAL output that is unique for your __samples__ (not found in blanks and QCs etc.)
- `intens_cols`: manual column indices of blanks, QCs, and samples.
- `info_cols`: manual column indices of the metabolite feature information (should be `1:32` for MS-DIAL version 4.70).
- `blankFC`: the fold-change increase compared to blank samples that features should have in order to be retained.
- `max_perc_mv`: samples with a greater percentage of missing values than this will be removed (e.g. `0.8`: samples with >80% missing values are removed).
- `missingPeaksFraction`: features with a greater percentage of missing peaks will be removed (e.g. `0.8`: features without peaks in >80% of samples are removed).
- `max_rsd`: set the maximum relative standard deviation of peaks that is acceptable within QC samples, and remove features that do not pass the test.
- `mv_imp_rowmax`: set the maximum acceptable percentage of missing data in any row.
- `mv_imp_colmax`: set the maximum acceptable percentage of missing data in any column.
- `mv_imp_method`: select the method you want to use for imputation of missing values. Most common options are k-nearest neighbours: `'knn'` and random forest: `'rf'`.

The default values for the arguments are shown here:

```r
pmp_preprocess <- function(pos_df, neg_df, metadata = NULL, samples_key = 'Sample', intens_cols = NULL, info_cols = NULL,
                           blankFC = 5, max_perc_mv = 0.8, missingPeaksFraction = 0.8, max_rsd = 25, 
                           mv_imp_rowmax = 0.7, mv_imp_colmax = 0.7, mv_imp_method = 'knn'){
                             # function code
                           }
```

#### Basic usage

At the simplest level, if you are satisfied with the default parameter values, you can implement function as shown below (an example for stool metabolomics).
It is still recommended that you manually set the columns containing the information and intensity values (of blanks, QCs, and samples) to ensure you get the results you expect.

```r
metab_pmp <- pmp_preprocess(pos_df = metab_stool_pos, neg_df = metab_stool_neg, samples_key = 'Stool',
                            intens_cols = c(33:54, 62:165), info_cols = 1:32)
```

You can then access the various list elements with standard list notation, e.g.:

```r
metab_stool_pmp$PCA_plot # view the PCA plot
metab_stool_pmp$glog_plot # confirm glog lambda value converged at the minima
metab_stool_pmp$filtering_dimensions # view how feature/samples numbers decreased with filtering
```

### Output

The function will return a list with five elements:

- `imputed_results`: a `SummarizedExperiment` object containing the PQN normalised imputed values.
- `glog_results`: a `SummarizedExperiment` object containining the PQN normalised, imputed, and glog-transformed values.
- `glog_plot`: a `ggplot2` plot showing the glog optimised lambda value to confirm that it converged at the minima.
- `PCA_plot`: a `ggplot2` PCA plot showing your samples by class (i.e. QC vs samples).
- `filtering_dimensions` a `data.frame` object showing the dimensions of your `SummarizedExperiment` object throughout filtering.
