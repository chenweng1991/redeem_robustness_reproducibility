library(ggplot2)
library(redeemR)
library(stringr)
library(dplyr)
library(gridExtra)
library(tidyr)
library(purrr)
library(Matrix)
library(Matrix.utils)
source("../utils/utils.r")

redeemv_final <- "/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/reproducibility_pub/data/redeemV_final/"
out <- "/lab/solexa_weissman/cweng/Projects/Collaborator/Caleb/github_redeem_plus/output/"

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
sample_folder <- args[2]


savename<-paste0(out,"/redeemr_bg.",name,".rds")


prepare_redeemr_bg<-function(redeemv_final, sample, reportpath){
    VariantsGTSummary<-redeemR.read(path=paste0(redeemv_final,sample),thr="S",Processed=F)
    redeemR<-Create_redeemR(VariantsGTSummary, maxctscut = 1,Cellcut = 1)
    redeemR<- add_raw_fragment(redeemR)
    saveRDS(redeemR,reportpath)
}

prepare_redeemr_bg(redeemv_final,sample_folder,savename)