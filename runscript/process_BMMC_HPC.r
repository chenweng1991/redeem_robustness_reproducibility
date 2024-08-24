library(redeemR)
library(dplyr)
source("../utils/utils.r")


args <- commandArgs(trailingOnly = TRUE)
sample1 <- args[1]
sample2 <- args[2]
name <- args[3]

WD <- "/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/reproducibility_pub/data/redeemV_final/"
output <- "../output/"
paths <- c(paste0(WD,sample1),paste0(WD,sample2))
savename<-paste0(output,"/combine_redeem_trim5.binom",name,".rds")
savename.GT<-paste0(output,"/combine_redeem_trim5.binom.GTSummary",name,".rds")

VariantsGTSummary<- redeemR.read.multiple.trim(paths=paths,thr="S", names=c("BMMC","HPC"), suffix= c("1", "2"),edge_trim=5)
# saveRDS(VariantsGTSummary,savename.GT)
print("GTSummary saved")


VariantsGTSummary <- readRDS(savename.GT)
redeemR<-Create_redeemR_model(VariantsGTSummary)
redeemR<- clean_redeem(redeemR,fdr = 0.05)
redeemR<-clean_redeem_removehot(redeemR)
#redeemR<- add_raw_fragment(redeemR,raw="RawGenotypes.Sensitive.StrandBalance")


savename<-paste0(output,"/combine_redeem_trim5.binom",name,".rds")
saveRDS(redeemR,savename)