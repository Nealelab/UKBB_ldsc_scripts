# /usr/bin/env Rscript

##################
#
# make manifest of UKB ldsc sumstat
# dropbox links 
#
##################

#### load data

# download of UKB GWAS file manifest
# https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=178908679
# manif <- readxl::read_excel("../reference/UKBB_GWAS_File_Manifest.xlsx",sheet = "Manifest 201807")
manif <- readxl::read_excel("../reference/UKBB GWAS Imputed v3 - File Manifest Release 20180731.xlsx",sheet = "Manifest 201807")
manif <- manif[grep("\\.gwas\\.imputed_v3\\.",manif$File),]

# ldsc results file
h2 <- read.delim("../results/round2_final/ukb31063_h2_all.02Oct2019.tsv.gz",header=T,stringsAsFactors=F,sep='\t')
h2 <- h2[,c("phenotype","sex","dilute","gwas_file","variable_type","source","description","isNotPrimary","confidence","h2_sig")]

# matching for biomarker gwas files with updated variant order 08/2021
h2$gwas_file[h2$source=="biomarkers" & !h2$dilute] <- stringr::str_replace(h2$gwas_file[h2$source=="biomarkers" & !h2$dilute], ".tsv.bgz", ".varorder.tsv.bgz")

# outer join on gwas file name
dat <- merge(h2, manif, by.x="gwas_file", by.y="File", all=TRUE)

# clean up matching fields with NAs from dilution, v2
# all.equal(dat$sex[!is.na(dat$sex) & !is.na(dat$Sex)],dat$Sex[!is.na(dat$sex) & !is.na(dat$Sex)])
# all.equal(dat$phenotype[!is.na(dat$sex) & !is.na(dat$Sex)],dat$`Phenotype Code`[!is.na(dat$sex) & !is.na(dat$Sex)])
dat$phenotype[is.na(dat$sex)] <- dat$`Phenotype Code`[is.na(dat$sex)]
dat$sex[is.na(dat$sex)] <- dat$Sex[is.na(dat$sex)]


# handle v2
# biomarkers duplicated in h2 but not manif due to dilution factor
dat$isV2 <- FALSE
dat$isV2[grepl(".v2.tsv",dat$gwas_file)] <- TRUE
biom_phens <- unique(dat$phenotype[dat$source=="biomarkers"])
dat$hasV2replacement <- FALSE
dat$hasV2replacement[duplicated(dat[,c("phenotype","sex")],fromLast=T) & !grepl(".v2.tsv",dat$gwas_file) & !(dat$phenotype %in% biom_phens)] <- TRUE
# table(grepl(".v2.tsv",dat$gwas_file),dat$hasV2replacement)
# all.equal(unique(dat$`Phenotype Description`[dat$isV2]),unique(dat$description[dat$isV2]))
dat$description[dat$hasV2replacement] <- dat$`Phenotype Description`[dat$hasV2replacement]
# all.equal(dat$`Phenotype Description`[dat$hasV2replacement],dat$description[dat$isV2][match(dat$phenotype[dat$hasV2replacement],dat$phenotype[dat$isV2])])
dat$source[dat$hasV2replacement] <- dat$source[dat$isV2][match(dat$phenotype[dat$hasV2replacement],dat$phenotype[dat$isV2])]
dat$variable_type[dat$hasV2replacement] <- dat$variable_type[dat$isV2][match(dat$phenotype[dat$hasV2replacement],dat$phenotype[dat$isV2])]

# fix icd showcase links
dat$`UK Biobank Data Showcase Link`[dat$source=="icd10"] <- "http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41202"
dat$`UK Biobank Data Showcase Link`[dat$source=="finngen"] <- "N/A"


# column cleanup
# dat <- dat[,-which(names(dat) %in% c("Sex","Phenotype Code","Phenotype Description"))]
names(dat)[names(dat)=="UK Biobank Data Showcase Link"] <- "showcase"

# confirm fields unique for mearge with dropbox files
# table(duplicated(dat[,c("phenotype","sex","isV2","dilute")]))


# ldsc dropbox file
box <- read.table("../reference/ukb31063_ldsc_sumstats_dropbox.txt",header=T,stringsAsFactors=F)

# parse info from filename
box$phenotype <- sapply(box$File, function(a) strsplit(a,'\\.')[[1]][1])
box$sex <- NA
box$sex[grep("\\.both_sexes\\.",box$File)] <- "both_sexes"
box$sex[grep("\\.female\\.",box$File)] <- "female"
box$sex[is.na(box$sex)] <- "male"
box$isV2 <- FALSE
box$isV2[grep("\\.v2\\.",box$File)] <- TRUE
box$dilute <- NA
box$dilute[grep("dilution_factor",box$File)] <- TRUE
box$dilute[!grepl("dilution_factor",box$File) & (box$phenotype %in% box$phenotype[!is.na(box$dilute) & box$dilute])] <- FALSE

# merge
# confirm fields unique for merge with dropbox files
# any(duplicated(dat[,c("phenotype","sex","isV2","dilute")]))
# any(duplicated(box[,c("phenotype","sex","isV2","dilute")]))
dat$mergeid <- paste0(dat$phenotype,"::",dat$sex,"::",dat$isV2,"::",dat$dilute)
box$mergeid <- paste0(box$phenotype,"::",box$sex,"::",box$isV2,"::",box$dilute)
any(duplicated(dat$mergeid))
any(duplicated(box$mergeid))
dat2 <- merge(dat,box,by="mergeid",all=TRUE)

# column cleanup
dat2$is_primary_gwas <- !dat2$isNotPrimary
names(dat2)[names(dat2)=="Dropbox"] <- "ldsc_sumstat_dropbox"
names(dat2)[names(dat2)=="AWS File"] <- "gwas_aws"
names(dat2)[names(dat2)=="wget command"] <- "gwas_wget"
names(dat2)[names(dat2)=="File"] <- "ldsc_sumstat_file"
dat2$ldsc_sumstat_wget <- paste0("wget ",dat2$ldsc_sumstat_dropbox," -O ",dat2$ldsc_sumstat_file)
names(dat2)[names(dat2)=="phenotype.x"] <- "phenotype"
names(dat2)[names(dat2)=="sex.x"] <- "sex"
names(dat2)[names(dat2)=="dilute.x"] <- "dilute"
names(dat2)[names(dat2)=="isV2.x"] <- "is_v2"
names(dat2)[names(dat2)=="hasV2replacement"] <- "has_v2_replacement"
dat2$is_primary_gwas[dat2$has_v2_replacement] <- FALSE
dat2$source[dat2$has_v2_replacement] <- "phesant"
names(dat2)[names(dat2)=="confidence"] <- "ldsc_confidence"
names(dat2)[names(dat2)=="h2_sig"] <- "ldsc_h2_significance"
dat2 <- dat2[,c("phenotype","description","showcase","source","sex","variable_type","is_v2","has_v2_replacement","dilute","ldsc_sumstat_file","ldsc_sumstat_dropbox","ldsc_sumstat_wget","gwas_file","gwas_aws","gwas_wget","is_primary_gwas","ldsc_confidence","ldsc_h2_significance")]

# save
con <- gzfile("../results/round2_final/ukb31063_ldsc_sumstat_manifest_aws_sep2022.tsv.gz","w")
write.table(dat2,file=con,col.names=T,row.names=F,sep='\t',quote=F)
close(con)

# eof