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



savename_trim5<-paste0(out,"/prepare_redeemr_trim5_v1.",name,".rds")
#savename_trim9<-paste0(out,"/prepare_redeemr_trim9_binom.",name,".rds")



prepare_redeemr_trim5_binom <- function(redeemv_final, sample, reportpath){
    VariantsGTSummary <- redeemR.read.trim(paste0(redeemv_final,sample), "S", Processed=F,
                                       "/VariantsGTSummary.S.trim5_binom.RDS",edge_trim=5)
    redeemR<-Create_redeemR_model(VariantsGTSummary)
    redeemR<- clean_redeem(redeemR,fdr = 0.05)
    redeemR<-clean_redeem_removehot(redeemR)
    redeemR<- add_raw_fragment(redeemR,raw="RawGenotypes.Sensitive.StrandBalance")
    redeem_qc_report<-run_redeem_qc(redeemR, redeemR@HomoVariants)
    report_trim4_binom<-list(redeemR=redeemR,redeem_qc_report=redeem_qc_report)
#    return(report_trim4_binom)
    saveRDS(report_trim4_binom,reportpath)
}

prepare_redeemr_trim9_binom <- function(redeemv_final, sample, reportpath){
    VariantsGTSummary <- redeemR.read.trim(paste0(redeemv_final,sample), "S", Processed=F,
                                       "/VariantsGTSummary.S.trim9_binom.RDS",edge_trim=9)
    redeemR<-Create_redeemR_model(VariantsGTSummary)
    redeemR<- clean_redeem(redeemR,fdr = 0.05)
    redeemR<-clean_redeem_removehot(redeemR)
    redeemR<- add_raw_fragment(redeemR,raw="RawGenotypes.Sensitive.StrandBalance")
    redeem_qc_report<-run_redeem_qc(redeemR, redeemR@HomoVariants)
    report_trim4_binom<-list(redeemR=redeemR,redeem_qc_report=redeem_qc_report)
#    return(report_trim4_binom)
    saveRDS(report_trim4_binom,reportpath)
}


## run the main function
prepare_redeemr_trim5_binom(redeemv_final,sample_folder,savename_trim5)
#prepare_redeemr_trim9_binom(redeemv_final,sample_folder,savename_trim9)
