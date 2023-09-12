#        +------------------------------------------------------+
#        |            matching SNPs for MarkerName              |
#        +------------------------------------------------------+
## import the PE gwas data and extract data via "Chromosome" variable one by one.
PE_gwas <- read.table("F:\\Doctorate\\aging_PE\\PE_GWAS_data\\preec.txt",
                      header = TRUE, sep = "\t")

chr22_gwas <- subset(PE_gwas, Chromosome == '22')
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
## install "BSgenome" and "SNPlocs.Hsapiens.dbSNP155.GRCh38" package through "BiocManager" package
#BiocManager::install("BSgenome")
#BiocManager::install("SNPlocs.Hsapiens.dbSNP11.GRCh11")
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
available.SNPs()

#get every Chromosome GRCh38 data vis "BSgenome" package
GRCh38_22 <- snpsBySeqname(SNPlocs.Hsapiens.dbSNP155.GRCh38, "22")
chr22_snps <- data.frame(GRCh38_22)

dim(chr22_snps)
write.csv(chr22_snps, "E:\\GRCh38\\chr22_match_snps.csv")
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# read the every Chromosome GRCh38 data into R.
chr22_snps <- read.table("F:\\GRCh38\\chr22_match_snps.csv",
                         header = TRUE, sep = ",")

## matching SNPs by "MarkerName"
chr22_snps$MarkerName <- paste(chr22_snps$seqnames,chr22_snps$pos,sep = ':')

chr22_gwas_snps <- merge(chr22_gwas,chr22_snps, by="MarkerName", all.x=TRUE)

dim(chr22_gwas_snps)

write.csv(chr22_gwas_snps, "F:\\Doctorate\\aging_PE\\PE_GWAS_data\\chr_gwas_snps\\chr16_gwas_snps.csv")
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
## after finishing the steps above for each Chromosome, All datasets should be set together.
# read every chromosome data matched with rsID.
library(vroom)
chr22_gwas_snps <- vroom("F:\\Doctorate\\aging_PE\\PE_GWAS_data\\chr_gwas_snps\\chr22_gwas_snps.csv",
                        delim = ",", col_names = TRUE)
# select the variables you needed.
library(dplyr)
pegwas22 <- select(chr22_gwas_snps, Chromosome,Position,MarkerName,RefSNP_id,Allele1,Allele2,
                  Freq1,FreqSE,MinFreq,MaxFreq,Effect,StdErr,P_value,Direction,
                  HetISq,HetChiSq,HetDf,HetPVal)
remove(chrX_gwas_snps)
# set data
PE_gwas <- rbind(pegwas1,pegwas2,pegwas3,pegwas4,pegwas5,pegwas6,pegwas7,
                 pegwas8,pegwas9,pegwas10,pegwas11,pegwas12,pegwas13,pegwas14,
                 pegwas15,pegwas16,pegwas17,pegwas18,pegwas19,pegwas20,
                 pegwas21,pegwas22,pegwas23,pegwasX)

write.table(PE_gwas,"F:\\Doctorate\\aging_PE\\PE_GWAS_data\\chr_gwas_snps\\PE_gwas.txt",
            sep = '\t', quote = F, row.names = F, col.names = T)