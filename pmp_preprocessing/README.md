# Peak Matrix Processing for raw MS-DIAL height matrices

The raw height matrices output by MS-DIAL will typically contain thousands (if not tens of thousands) of so-called "uninformative" and "irreproducible" features.
Features such as these could hinder downstream analyses, including statistical analysis, biomarker discovery, or pathway inference.

Therefore the common practice is to apply peak matrix validation and filtering procedures as described in Di Guida et al. 
([2016](https://doi.org/10.1007/s11306-016-1030-9)), Broadhurst et al. ([2018](https://doi.org/10.1007/s11306-018-1367-3)),
and Schiffman et al. ([2019](https://doi.org/10.1186/s12859-019-2871-9)).

Functions within the `pmp` (Peak Matrix Processing) package are designed to help users to prepare data for further statistical data analysis in a fast, 
easy to use and reproducible manner.

## Methods

An R function to handle the various steps is available within the [scripts folder](https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/tree/main/pmp_preprocessing/scripts), and will be explained in-depth below.

### R environment set-up

#### Load packages and display versions

Firstly, we need to load the required packages. We will also display their versions for future reference.

```r
# Get the R version
version$version.string

# Define a vector of required packages
pkgs <- c('here', 'data.table', 'tidyverse', 'kableExtra', 'ggplot2', 'pmp', 'SummarizedExperiment', 'S4Vectors')

# For each of the packages, load it and display the version number
for (pkg in pkgs) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  print(paste0(pkg, ': ', packageVersion(pkg)))
}

# Set the seed for generation of pseudo-random numbers
set.seed(2)
```

Example output:

```bash
[1] "R version 4.0.5 (2021-03-31)"
[1] "here: 1.0.1"
[1] "data.table: 1.14.0"
[1] "tidyverse: 1.3.1"
[1] "kableExtra: 1.3.4"
[1] "ggplot2: 3.3.5"
[1] "pmp: 1.2.1"
[1] "SummarizedExperiment: 1.20.0"
[1] "S4Vectors: 0.28.1"
```

### Data import

The data we will be importing is the height matrix data output from the MS-DIAL pipeline.
This output will be imported into a `SummarizedExperiment` object for pre-processing using the `pmp` package.

#### Metadata (optional at this stage)

_While this step is optional at this stage, the benefit of supplying metadata to the `SummarizedExperiment` object is that your metadata
will also be filtered by the pre-processing steps to match the remaining samples and features at the end._

Firstly, we will import and prepare the metadata component of our `SummarizedExperiment` object.
For simplicity and economy of code, we can specify the column types of our metadata using the `col_types` parameter in the `read_csv()` function.
These will usually be one of: `'c'` = character, `'f'` = factor, or `'n'` = numeric.

In this example, we will be working with LCMS data from a set of stool samples that were collected from different patients.

There are two metadata files: one that contains information about the samples themselves 
(e.g. sample name and ID, sampling site, and age at collection), and one that contains information about the patient 
(e.g. gender, disease type, other clinical metadata etc.). We will combine these two tables to produce a single metadata file for input
into the `SummarizedExperiment` object.

```r
# Import the patient metadata (the 'here' function makes defining file paths much easier)
metadata_patient <- read_csv(here::here('data', 'patient_metadata.csv'), col_types = 'ffffnnccfnc')

# Import the metabolomics sample metadata
metadata_metab_stool <- read_csv(here::here('data', 'metadata_metabolomics.csv'), col_types = 'ffffn') %>%
  filter(origin == 'stool') %>% # filter for the stool metabolomics metadata alone (not from other sites)
  left_join(metadata_patient, by = c('patient' = 'ID')) %>% # match patient metadata onto the sample metadata rows using a matched column
  column_to_rownames(var = 'sample') # choose the value that matches the MS-DIAL height matrix column names to be the rownames
```

#### Stool metabolomics data

We will load in both the positive and negative ionisation mode height matrices from the MS-DIAL pipeline.
We need to remove the MS/MS columns, as these were acquired in a different manner to the other columns.

Then, we can `rbind` the two tables together after denoting the original ionisation mode in the row names.
Because the QC data is incorporated within each feature (row), `pmp` can pre-process all of our data at once, and normalise the entire dataset.

```r
# Load the metabolomics height data
metab_stool_pos <- read_csv(here::here('data', 'stool_metabolites_height_positive.csv'), skip = 4) # skip the first 4 rows
metab_stool_neg <- read_csv(here::here('data', 'stool_metabolites_height_negative.csv'), skip = 4)

# Remove the MS/MS samples (not acquired in the same way)
metab_stool_pos <- metab_stool_pos[, !(names(metab_stool_pos) %in% c('MSMS_pos', 'MSMS_neg'))]
metab_stool_neg <- metab_stool_neg[, !(names(metab_stool_neg) %in% c('MSMS_pos', 'MSMS_neg'))]

# Separate into intensity and information data.frames for the SummarizedExperiment object
metab_stool_pos_counts <- as.matrix(metab_stool_pos[, c(33:54, 62:165)]) # you need to find the column numbers that contain Blanks, QCs, and samples
metab_stool_neg_counts <- as.matrix(metab_stool_neg[, c(33:54, 62:165)])

metab_stool_pos_info <- metab_stool_pos[, 1:32]
metab_stool_neg_info <- metab_stool_neg[, 1:32]

# Rename the data to indicate ionisation mode (using the MS-DIAL alignment ID)
stool_pos_rownames <- paste0(metab_stool_pos_info$`Alignment ID`, '_pos')
stool_neg_rownames <- paste0(metab_stool_neg_info$`Alignment ID`, '_neg')

rownames(metab_stool_pos_counts) <- stool_pos_rownames
rownames(metab_stool_pos_info) <- stool_pos_rownames
rownames(metab_stool_neg_counts) <- stool_neg_rownames
rownames(metab_stool_neg_info) <- stool_neg_rownames

# Merge the postive and negative ionisation modes using rbind
metab_stool_counts <- rbind(metab_stool_pos_counts, metab_stool_neg_counts)
metab_stool_info <- rbind(metab_stool_pos_info, metab_stool_neg_info)
```

Now that we have our intensity and feature information datasets, before we can import them into a `SummarizedExperiment` object,
we need to define class and group vectors so that `pmp` knows whether our variables are samples, QCs, or blanks.

The easiest way to achieve this is by getting the first two characters (a substring) of our column names,
and using these as indicators of the classes: `'Bl'` = blank, `'QC'` = QC etc.

```r
# Create class and group vectors
metab_stool_class <- substr(colnames(metab_stool_counts), start = 1, stop = 2)
metab_stool_group <- substr(colnames(metab_stool_counts), start = 1, stop = 2)
```

#### Import into a SummarizedExperiment object

Now that we have our individual components ready, we can check that the metadata matches our sample data,
and then import everything into a single, unified `SummarizedExperiment` object for pre-processing with the `pmp` package.

```r
# Reorder the metadata rows so that they match the column order of samples in the counts data.frame (ignoring blanks and QCs)
metadata_metab_stool <- metadata_metab_stool[colnames(metab_stool_counts)[23:126], ] # use only column containing sample intensities

# Check that the metadata matches the samples
identical(rownames(metadata_metab_stool), colnames(metab_stool_counts)[23:126]) # should return 'TRUE'

# Create the SummarizedExperiment (SE) object
metab_stool_SE <- SummarizedExperiment(assays = list(counts = metab_stool_counts),
                                       metadata = list(metadata_metab_stool),
                                       rowData = list(info = metab_stool_info),
                                       colData = DataFrame(class = metab_stool_class))
```

### Filtering and normalisation

#### Filtering on blank intensity

The first filtering steps we will carry out are to replace any so-called "missing" values (i.e. 0 values) with `NA` so
they will be compatible with downstream filtering steps.

After this, we can filter peaks based on the intensity values of our blanks.

```r
# Check the original number of features
dim(metab_stool_SE)

# Replace missing values with NA to be compatible with downstream filtering
assay(metab_stool_SE) <- replace(assay(metab_stool_SE), assay(metab_stool_SE) == 0, NA)

# Filter peaks and optionally samples based on blanks
metab_stool_filt <- filter_peaks_by_blanks(df = metab_stool_SE,
                                           fold_change = 5, # 5-fold change
                                           classes = metab_stool_SE$class,
                                           remove_samples = TRUE,
                                           remove_peaks = TRUE,
                                           blank_label = 'Bl') # from the class vector
                                           
# Check the number of features/samples remaining
dim(metab_stool_filt)
```

Example output (dimension checks):

```bash
[1] 81334   126
[1] 71700   115
```

#### Filtering on missing values and QC intensity variation

Next, we can perform a few filtering steps based on missing values and degree of variation in the QC samples.

- `filter_samples_by_mv`: removal of samples containing a user-defined maximum percentage of missing values [see documentation](https://rdrr.io/bioc/pmp/man/filter_samples_by_mv.html). E.g. setting the `max_perc_mv` argument to 0.8 will remove all samples with at least 80% missing values.
- `filter_peaks_by_fraction`: removal of peaks based upon the relative proportion of samples containing non-missing values [see documentation](https://rdrr.io/bioc/pmp/man/filter_peaks_by_fraction.html).
- `filter_peaks_by_rsd`: removal of peaks based upon relative standard deviation of intensity values for a given feature within specified QC samples [see documentation](https://rdrr.io/bioc/pmp/man/filter_peaks_by_rsd.html).

```r
# Filter samples based on the percentage of missing values
metab_stool_filt <- filter_samples_by_mv(df = metab_stool_filt,
                                         max_perc_mv = 0.8) # remove samples where >80% of values are missing 
                                         
# Check the number of features/samples
dim(metab_stool_filt)

# Filter peaks based on missing values across all samples
metab_stool_filt <- filter_peaks_by_fraction(df = metab_stool_filt,
                                             min_frac = 0.8, # features should have values for >80% of samples
                                             classes = metab_stool_filt$class,
                                             method = 'across')

# Check the number of features/samples
dim(metab_stool_filt)

# Filter peaks based on the percentage of variation in the QC samples
metab_stool_filt <- filter_peaks_by_rsd(df = metab_stool_filt,
                                        max_rsd = 25,
                                        classes = metab_stool_filt$class,
                                        qc_label = 'QC')

# Check the number of features/samples
dim(metab_stool_filt)
```

Example output (dimension checks):

```bash
[1] 71700   115
[1] 17315   115
[1] 11648   115
```

#### PQN normalisation and glog scaling

Finally we can normalise our data using probabilistic quotient normalisation (PQN) ([see documentation](https://rdrr.io/github/computational-metabolomics/pmp/man/pqn_normalisation.html)),
followed by missing value imputation (various method choices available; [see documentation](https://rdrr.io/github/computational-metabolomics/pmp/man/mv_imputation.html)), 
and then transform the data using a variance-stabling generalised logarithmic transformation ([see documentation](https://rdrr.io/github/computational-metabolomics/pmp/man/glog_transformation.html)).

The generalised logarithmic (glog) transformation stabilises the variance across low and high intensity mass spectral features.
The `glog_transformation` function uses QC samples to optimise the scaling factor `lambda`. 
Using the function `glog_plot_optimised_lambda` it’s possible to visualise if the optimisation of the given parameter has converged at the minima.

```r
# PQN data normalisation
metab_stool_norm <- pqn_normalisation(df = metab_stool_filt,
                                      classes = metab_stool_filt$class,
                                      qc_label = 'QC')
                                      
# Missing values imputation
metab_stool_imp <- mv_imputation(df = metab_stool_norm,
                                 rowmax = 0.7, # max % of missing data allowed in any row
                                 colmax = 0.7, # max % of missing data allowed in any column
                                 method = 'knn') # or rf, bcpa, sv, mn, md.
                                 
# Data scaling
metab_stool_glog <- glog_transformation(df = metab_stool_imp,
                                        classes = metab_stool_imp$class,
                                        qc_label = 'QC')

opt_lambda_stool <- processing_history(metab_stool_glog)$glog_transformation$lambda_opt

glog_plot_optimised_lambda(df = metab_stool_imp,
                           optimised_lambda = opt_lambda_stool,
                           classes = metab_stool_imp$class,
                           qc_label = 'QC')
```

Example plot:

<img src="https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/blob/main/pmp_preprocessing/assets/glog_optimised_lambda_plot.png" height="400">

We can also check the number of assigned metabolites that remain.

```r
# Number of assigned metabolites
table(rowData(metab_stool_glog)@listData[['info.Metabolite name']] != 'Unknown')
```

Example output (189 metabolites with assignments remain):

```bash
FALSE  TRUE 
11459   189 
```

### Principle component analysis

We can now visualise the output data using a principle component analysis plot (PCA) plot.

```r
# Perform the PCA and retrieve the explained variance values
PCA_stool <- prcomp(t(assay(metab_stool_glog)), center = TRUE)
varexp_stool <- c(summary(PCA_stool)$importance[2,1]*100,
                  summary(PCA_stool)$importance[2,2]*100)

# Create a dataset for plotting
data_PCA_stool <- cbind(data.frame(Samples = rownames(PCA_stool$x),
                                   PC1 = PCA_stool$x[,1],
                                   PC2 = PCA_stool$x[,2]),
                        class = metab_stool_glog$class)

# Plot results
(PCA_stool_plot <- ggplot(data_PCA_stool, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = class, color = factor(class))) +
    stat_ellipse(aes(fill = class), geom = 'polygon', type = 't', level = 0.9, alpha = 0.2) +
    labs(title = 'Stool Metabolomics',
         x = paste0('PC1 ', round(varexp_stool[1], 2), '%'),
         y = paste0('PC2 ', round(varexp_stool[2], 2), '%'))
)
```

Example plot:

<img src="https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/blob/main/pmp_preprocessing/assets/metab_stool_class_PCA.png" height="400">
