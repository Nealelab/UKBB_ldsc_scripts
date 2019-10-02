# /usr/bin/env Rscript

##################
#
# Goal: merge all h2 results to a single results file and
#       verify we have h2 for all current GWAS in manifest
#       with Ns matching the updated phenosummary files
#
##################

#####
# download of Neale UKB GWAS file manifest
#####

manif <- readxl::read_excel("../reference/UKBB_GWAS_File_Manifest.xlsx",sheet = "Manifest 201807")
gwas <- manif[grep("gwas",manif$File),]
rm(manif)

# remove entries for v1 gwases that have a v2
gwas <- gwas[!(duplicated(gwas[,c("Phenotype Code","Sex")],fromLast=T) & !grepl(".v2.tsv",gwas$File)),]


######
# load download of phenotype summary files
# see manifest
######

phensum <- read.delim("../reference/phenotypes.both_sexes.v2.tsv.bgz",header=T,stringsAsFactors=F,sep='\t')
phensum$sex <- "both_sexes"
phensum_male <- read.delim("../reference/phenotypes.male.v2.tsv.bgz",header=T,stringsAsFactors=F,sep='\t')
phensum_male$sex <- "male"
phensum_female <- read.delim("../reference/phenotypes.female.v2.tsv.bgz",header=T,stringsAsFactors=F,sep='\t')
phensum_female$sex <- "female"

# separate file for biomarkers
phenbio <- read.delim("../reference/biomarkers.both_sexes.tsv.bgz",header=T,stringsAsFactors=F,sep='\t')
phenbio$sex <- "both_sexes"
phenbio_male <- read.delim("../reference/biomarkers.male.tsv.bgz",header=T,stringsAsFactors=F,sep='\t')
phenbio_male$sex <- "male"
phenbio_female <- read.delim("../reference/biomarkers.female.tsv.bgz",header=T,stringsAsFactors=F,sep='\t')
phenbio_female$sex <- "female"


# combine files across sex conditions
phensum_all <- rbind(phensum,phensum_male,phensum_female,phenbio,phenbio_male,phenbio_female)
rm(phensum,phensum_male,phensum_female,phenbio_male,phenbio_female)
# keep phenbio for list of biomarker codes

######
# combine manifest files with meta-data
######

meta <- merge(gwas[,c("Phenotype Code","Phenotype Description","Sex","File")],
              phensum_all[,c("phenotype","description","variable_type","source","n_non_missing","n_missing","n_controls","n_cases","sex")],
              by.x=c("Phenotype Code","Sex"),
              by.y=c("phenotype","sex"))

# add entries for the biomarker alternative results with dilution factor (not in manifest)
meta$dilute <- NA
meta$dilute[meta$`Phenotype Code` %in% phenbio$phenotype] <- FALSE
metadil <- meta[meta$`Phenotype Code` %in% phenbio$phenotype,]
metadil$dilute <- TRUE
metadil$File <- sub("imputed_v3.","imputed_v3.dilution_factor.",metadil$File)
meta <- rbind(meta,metadil)

rm(gwas,phensum_all,phenbio,metadil)
# meta[trimws(meta$description) != trimws(meta$`Phenotype Description`),c("description","Phenotype Description")]

######
# load download of phesant h2 results from google bucket
######

dat <- read.delim("../results/round2_raw/ukbb31063.both_sexes.h2part_results.v2.phesant.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat$sex <- "both_sexes"
dat_male <- read.delim("../results/round2_raw/ukbb31063.male.h2part_results.v2.phesant.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_male$sex <- "male"
dat_female <- read.delim("../results/round2_raw/ukbb31063.female.h2part_results.v2.phesant.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_female$sex <- "female"
dat_phes_all <- rbind(dat,dat_male,dat_female)
dat_phes_all$dilute <- NA
rm(dat,dat_male,dat_female)


######
# load download of icd10, finngen, covariate, biomarker h2 results from google bucket
######

dat_icd <- read.delim("../results/round2_raw/ukbb31063.both_sexes.h2part_results.icd10.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_icd$sex <- "both_sexes"
dat_icd_female <- read.delim("../results/round2_raw/ukbb31063.female.h2part_results.icd10.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_icd_female$sex <- "female"
dat_icd_male <- read.delim("../results/round2_raw/ukbb31063.male.h2part_results.icd10.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_icd_male$sex <- "male"
dat_icd_all <- rbind(dat_icd,dat_icd_female,dat_icd_male)
dat_icd_all$dilute <- NA
rm(dat_icd,dat_icd_male,dat_icd_female)

dat_fin <- read.delim("../results/round2_raw/ukbb31063.both_sexes.h2part_results.finngen.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_fin$sex <- "both_sexes"
dat_fin_female <- read.delim("../results/round2_raw/ukbb31063.female.h2part_results.finngen.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_fin_female$sex <- "female"
dat_fin_male <- read.delim("../results/round2_raw/ukbb31063.male.h2part_results.finngen.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_fin_male$sex <- "male"
dat_fin_all <- rbind(dat_fin,dat_fin_female,dat_fin_male)
dat_fin_all$dilute <- NA
rm(dat_fin,dat_fin_male,dat_fin_female)

dat_cov <- read.delim("../results/round2_raw/ukbb31063.both_sexes.h2part_results.covariate.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_cov$sex <- "both_sexes"
dat_cov_female <- read.delim("../results/round2_raw/ukbb31063.female.h2part_results.covariate.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_cov_female$sex <- "female"
dat_cov_male <- read.delim("../results/round2_raw/ukbb31063.male.h2part_results.covariate.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_cov_male$sex <- "male"
dat_cov_all <- rbind(dat_cov,dat_cov_female,dat_cov_male)
dat_cov_all$dilute <- NA
rm(dat_cov,dat_cov_male,dat_cov_female)

dat_bio <- read.delim("../results/round2_raw/ukbb31063.both_sexes.h2part_results.biomarkers.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_bio$sex <- "both_sexes"
dat_bio_female <- read.delim("../results/round2_raw/ukbb31063.female.h2part_results.biomarkers.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_bio_female$sex <- "female"
dat_bio_male <- read.delim("../results/round2_raw/ukbb31063.male.h2part_results.biomarkers.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_bio_male$sex <- "male"
dat_bio_all <- rbind(dat_bio,dat_bio_female,dat_bio_male)
dat_bio_all$dilute <- FALSE
rm(dat_bio,dat_bio_male,dat_bio_female)

dat_biod <- read.delim("../results/round2_raw/ukbb31063.both_sexes.h2part_results.biomarkers.dilution.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_biod$sex <- "both_sexes"
dat_biod_female <- read.delim("../results/round2_raw/ukbb31063.female.h2part_results.biomarkers.dilution.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_biod_female$sex <- "female"
dat_biod_male <- read.delim("../results/round2_raw/ukbb31063.male.h2part_results.biomarkers.dilution.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_biod_male$sex <- "male"
dat_biod_all <- rbind(dat_biod,dat_biod_female,dat_biod_male)
dat_biod_all$dilute <- TRUE
rm(dat_biod,dat_biod_male,dat_biod_female)

# merge
h2 <- rbind(dat_phes_all,dat_icd_all,dat_fin_all,dat_cov_all,dat_bio_all,dat_biod_all)
rm(dat_phes_all,dat_icd_all,dat_fin_all,dat_cov_all,dat_bio_all,dat_biod_all)


########
# align with meta-data
########

dat <- merge(meta,h2,by.x=c("Phenotype Code","Sex","dilute"),by.y=c("phenotype","sex","dilute"),all=T)

### some sanity checks
# all.equal(dat$n,dat$n_non_missing)
# all.equal(dat$n_cases.x,dat$n_cases.y)
# all.equal(dat$n_controls.x,dat$n_controls.y)
# all.equal(dat$source.x,dat$source.y) # false, in the x version is NA for "covariate" phenotypes
# dat[dat$description.x != dat$description.y, c("Phenotype Description","description.x","description.y")] # finngen NAs, and minor naming
# dat[dat$description.y != dat$`Phenotype Description`, c("Phenotype Description","description.x","description.y")] # additional vitamin/mineral
# dat[trimws(dat$description.y) != trimws(dat$`Phenotype Description`), c("Phenotype Description","description.x","description.y")] # fixes vitamin


#######
# Update names of finngen phenos
#######

# file from Duncan Palmer with updated descriptions
fin_name <- read.table("../reference/ukb_finngen_names.tsv",sep='\t',comment.char = "", header=T, stringsAsFactors = F)

fin_idx <- match(dat$`Phenotype Code`,fin_name$NAME)
dat$description.y[!is.na(fin_idx)] <- fin_name$LONGNAME[fin_idx[!is.na(fin_idx)]]

# any(dat$description.y=="nan")


#######
# clean and save
#######

dat <- dat[,!(names(dat) %in% c("source.x","n_cases.x","n_controls.x","n_non_missing","description.x","Phenotype Description"))] # duplicated
names(dat)[names(dat)=="n_cases.y"] <- "n_cases"
names(dat)[names(dat)=="n_controls.y"] <- "n_controls"
names(dat)[names(dat)=="description.y"] <- "description"
names(dat)[names(dat)=="source.y"] <- "source"
names(dat)[names(dat)=="Phenotype Code"] <- "phenotype"
# dat$description[dat$description=="nan"] <- dat$phenotype[dat$description=="nan"]
dat$description <- trimws(dat$description)
names(dat)[names(dat)=="Sex"] <- "sex"
names(dat)[names(dat)=="File"] <- "gwas_file"

# make biomarker name consistent between raw/irnt verison
dat$description[endsWith(dat$description,"(quantile)")] <- dat$description[match(sub("_irnt","_raw",dat$phenotype[endsWith(dat$description,"(quantile)")]),dat$phenotype)]


# fix p=0
require("Rmpfr")
for(jj in grep("_p$",names(dat))){
  if(any(!is.na(dat[,jj]) & dat[,jj]==0)){
    jjz <- grep(gsub("_p$","_z",names(dat)[jj]),names(dat))
    print(paste("replacing p=0 in",names(dat)[jj],"using z score",names(dat)[jjz],sep=" "))
    zero_p <- which(!is.na(dat[,jj]) & dat[,jj]==0)
    dat[zero_p,jj] <- format(mpfr(2,64)*pnorm(mpfr(dat[[jjz]][zero_p],64),lower=F), max.digits=15, scientific=T)
    
    # alt method that would extend even further
    # and could add some nchar check here to get fixed width on the whole thing
    # lp <- (pnorm(mpfr(dat[[jjz]][zero_p]),lower=F, log.p=T)/log(10)) + log10(2)
    # flp <- floor(lp)
    # dat[zero_p,jj] <- paste0(format(10^(lp-flp),max.digits=6,scientific=F),"e",sprintf("%03d",flp))
    
    }
}
head(dat$intercept_p[dat$intercept_z > 40])

con <- gzfile("../results/round2_final/ukb31063.h2.baseline1-1.gwas_v2.full_results.02Oct2019.tsv.gz","w")
write.table(dat,file=con,sep='\t',col.names=T,row.names=F,quote=F)
close(con)

