#! /usr/bin/env bash

####
#
# Download raw results and reference files for
# the UKB h2 analyses
#
####

# Neale lab GWAS manifest
curl "https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/export?format=xlsx" --output "../reference/UKBB_GWAS_File_Manifest.xlsx"


# phenotype files
# wget https://www.dropbox.com/s/d4mlq9ly93yhjyt/phenotypes.both_sexes.tsv.bgz -O ../reference/phenotypes.both_sexes.tsv.gz
# wget https://www.dropbox.com/s/r7idtoulxfpyjss/phenotypes.female.tsv.bgz -O ../reference/phenotypes.female.tsv.gz
# wget https://www.dropbox.com/s/mywldevz4nsla2r/phenotypes.male.tsv.bgz -O ../reference/phenotypes.male.tsv.gz
wget https://www.dropbox.com/s/j1h9e0lblmukiko/phenotypes.both_sexes.v2.tsv.bgz?dl=0 -O ../reference/phenotypes.both_sexes.v2.tsv.bgz
wget https://www.dropbox.com/s/0a87i4y049vje12/phenotypes.female.v2.tsv.bgz?dl=0 -O ../reference/phenotypes.female.v2.tsv.bgz
wget https://www.dropbox.com/s/k8h7j85awav0lrt/phenotypes.male.v2.tsv.bgz?dl=0 -O ../reference/phenotypes.male.v2.tsv.bgz


# additional biomarker phenotypes
wget https://www.dropbox.com/s/vuzx5bkv0dfwtn7/biomarkers.both_sexes.tsv.bgz?dl=0 -O ../reference/biomarkers.both_sexes.tsv.bgz
wget https://www.dropbox.com/s/f549pns6tb3su8r/biomarkers.female.tsv.bgz?dl=0 -O ../reference/biomarkers.female.tsv.bgz
wget https://www.dropbox.com/s/vmq6ykz0ld74r2n/biomarkers.male.tsv.bgz?dl=0 -O ../reference/biomarkers.male.tsv.bgz


# names for finngen codes
# ukb_finngen_names.tsv (internal file currently)


# h2 results
gsutil cp gs://nbaya/h2part/updates/ukbb31063.*.h2part_results.v2.phesant.tsv.gz ../results/round2_raw/

for sex in both_sexes female male; do
	gsutil cp gs://nbaya/h2part/ukbb31063.${sex}.h2part_results.icd10.tsv.gz ../results/round2_raw/
	gsutil cp gs://nbaya/h2part/ukbb31063.${sex}.h2part_results.finngen.tsv.gz ../results/round2_raw/
	gsutil cp gs://ukb31063-ldsc-results/h2part_misc/ukbb31063.${sex}.h2part_results.covariate.batch_1.tsv.gz ../results/round2_raw/
	gsutil cp gs://ukb31063-ldsc-results/h2part_misc/ukbb31063.${sex}.h2part_results.biomarkers.batch_1.tsv.gz ../results/round2_raw/
	gsutil cp gs://ukb31063-ldsc-results/h2part_misc/ukbb31063.${sex}.h2part_results.biomarkers.dilution.batch_1.tsv.gz ../results/round2_raw/
done


# phenotypic correlations
scp rwalters@login.broadinstitute.org:/stanley/robinson/ccarey/UKBB/h2_corrmat/corrmat/*resid_corrmat.csv ../reference/

# pairwise complete sample sizes
scp rwalters@login.broadinstitute.org:/stanley/robinson/ccarey/UKBB/h2_corrmat/corrmat/*pairwise_complete_ns.csv ../reference/

# UKB phenotype codings
wget https://github.com/astheeggeggs/PHESANT/raw/master/variable-info/Data_Dictionary_Showcase.csv -O ../reference/Data_Dictionary_Showcase.csv

# identified suboptimal ordinal codes
# ukb_ord_codings_warn.txt (manually added)

# eof
