
##################
#
# Goal: verify we have h2 for all current GWAS in manifest
#       with Ns matching the updated phenosummary files
#
##################


#####
# download of UKB GWAS file manifest
# https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=178908679
#####

manif <- readxl::read_excel("UKBB GWAS Imputed v3 - File Manifest Release 20180731.xlsx",sheet = "Manifest 201807")
gwas <- manif[grep("gwas",manif$File),]
rm(manif)

# remove entries for v1 gwases that have a v2
gwas <- gwas[!(duplicated(gwas[,c("Phenotype Code","Sex")],fromLast=T) & !grepl(".v2.tsv",gwas$File)),]


######
# load download of phenotype summary files
# see manifest
######

phensum <- read.delim("phenotypes.both_sexes.v2.tsv.bgz",header=T,stringsAsFactors=F,sep='\t')
phensum$sex <- "both_sexes"
phensum_male <- read.delim("phenotypes.male.v2.tsv.bgz",header=T,stringsAsFactors=F,sep='\t')
phensum_male$sex <- "male"
phensum_female <- read.delim("phenotypes.female.v2.tsv.bgz",header=T,stringsAsFactors=F,sep='\t')
phensum_female$sex <- "female"

# combine files across sex conditions
phensum_all <- rbind(phensum,phensum_male,phensum_female)
rm(phensum,phensum_male,phensum_female)


######
# combine manifest files with meta-data
######

meta <- merge(gwas[,c("Phenotype Code","Phenotype Description","Sex","File")],
              phensum_all[,c("phenotype","description","variable_type","source","n_non_missing","n_missing","n_controls","n_cases","sex")],
              by.x=c("Phenotype Code","Sex"),
              by.y=c("phenotype","sex"))

rm(gwas,phensum_all)
# meta[trimws(meta$description) != trimws(meta$`Phenotype Description`),c("description","Phenotype Description")]

######
# load download of update phesant h2 results from google bucket
# gsutil cp gs://nbaya/h2part/updates/ukbb31063.*.h2part_results.v2.phesant.tsv.gz ./
######

dat <- read.delim("ukbb31063.both_sexes.h2part_results.v2.phesant.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat$sex <- "both_sexes"
dat_male <- read.delim("ukbb31063.male.h2part_results.v2.phesant.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_male$sex <- "male"
dat_female <- read.delim("ukbb31063.female.h2part_results.v2.phesant.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_female$sex <- "female"
dat_phes_all <- rbind(dat,dat_male,dat_female)
rm(dat,dat_male,dat_female)


######
# load download of icd10, finngen, covariate h2 results from google bucket (no updates)
# gsutil cp gs://nbaya/h2part/ukbb31063.{both_sexes,female,male}.h2part_results.{icd10,finngen}.tsv.gz ./
# gsutil cp gs://ukb31063-ldsc-results/h2part_misc/ukbb31063.{both_sexes,female,male}.h2part_results.covariate.batch_1.tsv.gz ./
######

dat_icd <- read.delim("ukbb31063.both_sexes.h2part_results.icd10.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_icd$sex <- "both_sexes"
dat_icd_female <- read.delim("ukbb31063.female.h2part_results.icd10.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_icd_female$sex <- "female"
dat_icd_male <- read.delim("ukbb31063.male.h2part_results.icd10.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_icd_male$sex <- "male"
dat_icd_all <- rbind(dat_icd,dat_icd_female,dat_icd_male)
rm(dat_icd,dat_icd_male,dat_icd_female)

dat_fin <- read.delim("ukbb31063.both_sexes.h2part_results.finngen.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_fin$sex <- "both_sexes"
dat_fin_female <- read.delim("ukbb31063.female.h2part_results.finngen.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_fin_female$sex <- "female"
dat_fin_male <- read.delim("ukbb31063.male.h2part_results.finngen.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_fin_male$sex <- "male"
dat_fin_all <- rbind(dat_fin,dat_fin_female,dat_fin_male)
rm(dat_fin,dat_fin_male,dat_fin_female)

dat_cov <- read.delim("ukbb31063.both_sexes.h2part_results.covariate.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_cov$sex <- "both_sexes"
dat_cov_female <- read.delim("ukbb31063.female.h2part_results.covariate.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_cov_female$sex <- "female"
dat_cov_male <- read.delim("ukbb31063.male.h2part_results.covariate.batch_1.tsv.gz",header=T,stringsAsFactors=F,sep='\t',quote="")
dat_cov_male$sex <- "male"
dat_cov_all <- rbind(dat_cov,dat_cov_female,dat_cov_male)
rm(dat_cov,dat_cov_male,dat_cov_female)

# merge
h2 <- rbind(dat_phes_all,dat_icd_all,dat_fin_all,dat_cov_all)
rm(dat_phes_all,dat_icd_all,dat_fin_all,dat_cov_all)

########
# align with meta-data
########

dat <- merge(meta,h2,by.x=c("Phenotype Code","Sex"),by.y=c("phenotype","sex"),all=T)

### some sanity checks
# all.equal(dat$n,dat$n_non_missing)
# all.equal(dat$n_cases.x,dat$n_cases.y)
# all.equal(dat$n_controls.x,dat$n_controls.y)
# all.equal(dat$source.x,dat$source.y) $ false, is the x version is NA for "covariate" phenotypes
# dat[dat$description.x != dat$description.y, c("Phenotype Description","description.x","description.y")] # finngen NAs, and minor naming
# dat[dat$description.y != dat$`Phenotype Description`, c("Phenotype Description","description.x","description.y")] # additional vitamin/mineral
# dat[trimws(dat$description.y) != trimws(dat$`Phenotype Description`), c("Phenotype Description","description.x","description.y")] # fixes vitamin


#######
# Update names of finngen phenos
#######

# file from Duncan with updated descriptions
fin_name <- read.table("ukb_finngen_names.tsv",sep='\t',comment.char = "", header=T, stringsAsFactors = F)

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


# fix p=0
require("Rmpfr")
for(jj in grep("_p$",names(dat))){
  if(any(!is.na(dat[,jj]) & dat[,jj]==0)){
    jjz <- grep(gsub("_p$","_z",names(dat)[jj]),names(dat))
    print(paste("replacing p=0 in",names(dat)[jj],"using z score",names(dat)[jjz],sep=" "))
    zero_p <- which(!is.na(dat[,jj]) & dat[,jj]==0)
#    dat[[jj]] <- mpfr(dat[[jj]],64)
#    dat[[jj]][zero_p] <- mpfr(0.5,64)*erfc(mpfr(dat[[jjz]][zero_p],64)/sqrt(mpfr(2,64)))
    dat[zero_p,jj] <- format(mpfr(0.5,64)*erfc(mpfr(dat[[jjz]][zero_p],64)/sqrt(mpfr(2,64))), max.digits=15, scientific=T)
    }
}
head(dat$intercept_p[dat$intercept_z > 40])

con <- gzfile("ukb31063.h2.baseline1-1.gwas_v2.full_results.tsv.gz","w")
write.table(dat,file=con,sep='\t',col.names=T,row.names=F,quote=F)
close(con)

