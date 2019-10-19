# UKBB ldsc scripts

Scripts used for running LD score regression ([LDSC](https://github.com/bulik/ldsc)) to estimate SNP heritability in UK Biobank. See the [Neale Lab UK Biobank GWAS](http://www.nealelab.is/uk-biobank), the [heritability analysis post](http://www.nealelab.is/blog/2017/9/15/heritability-of-2000-traits-and-disorders-in-the-uk-biobank) and the [results site](https://nealelab.github.io/UKBB_ldsc/).

These scripts are mostly intended for documentation purposes, and may take some work to get running or adapt to other applications. The code used to generate the [results site](https://nealelab.github.io/UKBB_ldsc/) is maintained separately in the [UKBB_ldsc repo](github.com/nealelab/UKBB_ldsc).

# Table of Contents

* [Usage](#usage)
  * [Requirements](#requirements)
  * [Running ldsc](#running-ldsc)
    * [Implementation Note](#implementation-note)
    * [Settings](#settings)
    * [Submission script](#submission-script)
  * [Aggregating results](#aggregating-results)
  * [Processing results](#processing-results)
  * [Preparing for results site](#preparing-for-results-site)
* [Results](#results)
* [Round 1 archive](#round-1-archive)
* [Credits/Contact](#credits)

TOC originally created with [gh-md-toc](https://github.com/ekalinin/github-markdown-toc.go)


# Usage

## Requirements

Broadly, the core scripts here depend on:

* [Google Cloud Platform](https://cloud.google.com/), for storage and compute
* [cloudtools](https://github.com/Nealelab/cloudtools), for interaction with google dataproc
 * Run prior to the introduction of [hailctl](https://hail.is/docs/0.2/hail_on_the_cloud.html), which is now preferred
* [Hail](https://hail.is/), for cloud compute support
* [MTAG](https://github.com/omeed-maghzian/mtag), for it's modified implementation of [LDSC](https://github.com/bulik/ldsc)

See the respective scripts for more info.


## Setup 

The scripts here assume [LDSC-formatted sumstats files](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format) have already been created. For the current UKB GWAS, that process is handled by code in the (UK_Biobank_GWAS repo)[https://github.com/Nealelab/UK_Biobank_GWAS], specifically `./0.1/27.export_ldsc_sumstats.py` and `./0.2/export_ldsc_sumstats_biomarkers.py` (for Hail versions 0.1 and 0.2, respectively). The creation of a Hail table listing HapMap3 sites passing QC in UK Biobank to be extracted is also now covered in that repo by `./0.1/26.create_ldsc_hm3_keytable.py` and `./0.2/create_ldsc_hm3_table.py`. Previous code for that process in the Round 1 analysis is in the `./round1/` archive folder of this repo.

Copies of the LDSC-formatted sumstats files used for the current analysis are [available for download](https://nealelab.github.io/UKBB_ldsc/douwloads.html)


## Running ldsc

* Job scripts: `ldsc_h2part_parallel_batch_v2*.py`
* Submission script: `ldsc_h2part_parallel_batch_submit.sh`

These script run [ldsc](https://github.com/bulik/ldsc) to estimate SNP-heritability for a batch of phenotypes using [MTAG](https://github.com/omeed-maghzian/mtag). We use the version of ldsc provided with MTAG to take advantage of it's interface for calling ldsc from within python rather than via the command line. 

LD score regression is run using the `baselineLD_v1.1` partitioned model with the precomputed LD scores from 1000 Genomes Phase 3 European populations available (here)[https://data.broadinstitute.org/alkesgroup/LDSCORE/] and described by [Gazal et al. (2017)](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3954.html). The full set of results for the annotations (e.g. corresponding to the `--print-coefficients` output file) are converted to a flat row of results per phenotype. Previous code for running univariate LDSR is in the `./round1/` directory.

For binary traits heritability results on the liability scale are also provided. Conversion to liability-scale h2 for dichotomous traits is done assuming that the population prevelence is equal to the prevelence in the UKBB GWAS sample. 

Versions of the job script (e.g. `biomarkers`, `misc`) differ slightly due to slight differences in file naming conventions and other settings across batch of the UK Biobank GWAS and LDSC-formatted file export, but all function the same way.

### Settings

Using this script involves both setting job information in the script and passing arguments for parallelization via the command line. The following setting are set by editing the the python script (`ldsc_h2part_parallel_batch_v2.py`):

```
wd = '/home/mtag/' # root working directory
ld_ref_panel = '/home/mtag/mtag-master/ld_ref_panel/baselineLD_v1.1/baselineLD.' # local path
ld_w_panel = '/home/mtag/mtag-master/ld_ref_panel/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.' # local path
ld_frq_panel = '/home/mtag/mtag-master/ld_ref_panel/1000G_Phase3_frq/1000G.EUR.QC.' # local path
phen_summary = 'gs://ukbb-gwas-imputed-v3-results/export2/phenotypes.{}.tsv.bgz' # google cloud bucket file
ss_bucket = 'gs://ukb31063-mega-gwas/ldsc-export/sumstats-files' # google cloud bucket folder
out_bucket = 'gs://ukb31063-ldsc-results/h2part_misc/' # ouput google bucket location
num_proc = 6 # number of processes to run
```

The three LD panel arguments are paths and file stems to the LD score reference files locally on the GCP node, such that this string can be supplied directly to the corresponding argument of [ldsc](https://github.com/bulik/ldsc). If the files are missing the script will download the above-named default files, but the code will require modification if you want it to download some other LD reference.

The structure of the `phen_summary` file is assumed to match the phenotype files provided with the [Neale UKB GWAS](http://www.nealelab.is/uk-biobank). Downloads avaialble [here](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU).

From the command line, this script expects:

* `--sex-group STRING`, whether the the GWAS is for `female`,`male`, or `both_sexes` (because it affects our ldsc file names)
* `--phsource STRING`, e.g. `phesant` or `icd10`, to extract for the `phen_summary` file
* `--numphens INT`, the total number of phenotypes in this category, for the purpose of parallelization
* `--parsplit INT`, the number of parallel batches to split phenotypes into
* `--paridx INT`, the index of the current batch amoung the `parsplit` batches

It's assumed that there are `numphens` phenotypes in the `phen_summary` file with source `phsource` and they they all have ldsc sumstats results files in the specified bucket. The mapping between the phenotypes names and the filenames for the sumstats files is hardcoded in the script (see `ss_name`).

Within an instance of `ldsc_h2part_parallel_batch.py`, it will estimate h2 for `num_proc` phenotypes at a time in parallel, and continue running until it's share of the `num_phens` phenotypes is finished.


### Submission script

 Script: `ldsc_h2*_parallel_batch_submit.sh`

**WARNING: THIS SCRIPT WILL CREATE NEW CLUSTERS BUT WILL NOT SHUT THEM DOWN**

This is a simple bash script to loop job submission of `ldsc_h2part_parallel_batch_v2*py` for each batch of phenotypes, assuming parallelization to `$maxi` batches. A new cluster is spun up and a job submitted to each cluster using [cloudtools](https://github.com/Nealelab/cloudtools). If using this script, please monitor the jobs and **shut down each cluster when it completes**.  
Although we document here the submission script we used based on [cloudtools](https://github.com/Nealelab/cloudtools), the same structure can be applied using [hailctl](https://hail.is/docs/0.2/hail_on_the_cloud.html) now provided with [Hail](https://hail.is/). We recommend that approach for future applications.

### Implementation Note

This is almost certainly not the ideal way to structure this analysis. Making a human manage the splitting/batching here somewhat defeats the purpose of having flexible cloud compute. We're only doing it this way currently for expediency.

With the current settings in the submission script running univariate h2 for ~1500 traits takes a little over 10 node hours (~160 CPU hours) split over 10 minimal `n1-highcpu-16` VMs, each running 8 traits in parallel at a time. Paritioned heritbaility is slower, taking just under 20 hours on 10 `n1-standard-16` VMs running 6 traits at a time for the same ~1500 traits. Attempts to scale this up to more traits on a machine (with 32- or 64-core VMs) have seen poor performance, likely due to being I/O bound for reading reference and sumstat files. 


## Aggregating results

Scripts: `fetch_files.sh`, `ukbb_h2_results_merge.R`

`fetch_files.sh` downloads the completed LD score regression results and other reference files used for processing the results. Note that the output locations go to other folders (based on if these scripts are a subfolder within the [site repo structure](https://github.com/Nealelab/UKBB_ldsc)). Not all of the listed files are currently publicly available, most of the primary ones are.

`ukbb_h2_results_merge.R` merges the different batches of LD score regression results. Descriptive information is also incorporated from the appropriate phenosummary files and the GWAS file manifest. 

These scripts replace the previous `./round1/agg_ldsc.sh`, though as is evident from the file names in these scripts we did some aggregation prior to writing these cleaner scripts. These should however be generalizeable going forward.


## Processing results

Scripts: `_ukbb_h2_select_topline.Rmd`, `_ukbb_h2_process_confidence.Rmd`, `_ukbb_h2_process_significance.Rmd`

From the LD score results for the full set of round 2 GWAS, we perform several steps of post-processing documented in these R markdown scripts. The full contents of these steps are included in the [methods section of the results site](https://nealelab.github.io/UKBB_ldsc/details.html). Briefly, these steps are as follows:

* `_ukbb_h2_select_topline.Rmd`: Select a primary GWAS for each unqiue phenotype in the Neale Lab GWAS. This includes choices of raw vs. rank-normalized versions of continuous phenotypes, whether to add a covariate for the estimated dilution fraction for biomarker phenotypes, and whether some phenotypes need sex-specific results. The output of this process is available [here](https://nealelab.github.io/UKBB_ldsc/select_topline.html).

* `_ukbb_h2_process_confidence.Rmd`: Define confidence ratings for the LDSR results. This largely involves consider power and biases as a function of sample size, but also considers the behavior of the standard errors and possible issues with phenotype codings. The output of this process is available [here](https://nealelab.github.io/UKBB_ldsc/confidence.html).

* `_ukbb_h2_process_significance.Rmd`: Identify statistically significant SNP heritability results. This process focuses on evaluating the multiple testing burden for the thousands of available UK Biobank GWAS. The output of this process is available [here](https://nealelab.github.io/UKBB_ldsc/significance.html).

Note: These Rmd files are designed to compile as part of the [results site code](https://github.com/Nealelab/UKBB_ldsc) rather than as stand-alone files here. As such they rely on the `_code_highlight_fix.Rmd` and `_toc_fix.Rmd` fragments, which are available in `./site_source/rmd_source/` of the [UKBB_ldsc repo](https://github.com/Nealelab/UKBB_ldsc), and are best compiled using the CSS files from that repo.


## Preparing for results site

Script: `build_ldsc_sumstat_manifest.R`

This script builds a manifest of the downloads links for the LDSR sumstat files that are [available for download](https://nealelab.github.io/UKBB_ldsc/downloads.html). This includes aligning those links with the GWAS manifest and providing additional fields for convenient filtering and programmatic access to the files (e.g. `wget` commands). Dropbox links were compiled outside of this script with the assistance of [Dropbox Uploader](https://github.com/andreafabrizi/Dropbox-Uploader).


# Results

* [Results browser](https://nealelab.github.io/UKBB_ldsc/h2_browser.html)
* [Full downloads](https://nealelab.github.io/UKBB_ldsc/downloads.html)

For more about the underlying GWAS, see [nealelab.is/uk-biobank](http://www.nealelab.is/uk-biobank)

# Round 1 archive

Scripts for the LD Score regression results for the previous round of UK Biobank GWAS are available in the `./round1/` folder. The corresponding results are [archived on the results site](https://nealelab.github.io/UKBB_ldsc/round1_h2_browser.html) and remain [available for download](https://nealelab.github.io/UKBB_ldsc/docs/downloads.html#previous_snp_heritability_results).


# Credits

* Lead: Raymond Walters
* LDSR analysis support: Nikolas Baya, Katherine Tashman, Danfeng Chen, Liam Abbott
* Phenotype analysis support: Caitlin Carey, Duncan Palmer
* Supervision: Benjamin Neale
* UK Biobank GWAS Core Team: Liam Abbott, Sam Bryant, Claire Churchhouse, Andrea Ganna, Daniel Howrigan, Duncan Palmer, Ben Neale, Raymond Walters, Caitlin Carey, The Hail team, Benjamin Neale

For full credits see the [results site](https://nealelab.github.io/UKBB_ldsc/credits.html).

For any questions or feedback, please contact us at: **nealelab.ukb@gmail.com**
