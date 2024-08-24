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


savename_1<-paste0(out,"/prepare_redeemr_v1.",name,".rds")
savename_2<-paste0(out,"/prepare_redeemr_trim4_v1.",name,".rds")
savename_3<-paste0(out,"/prepare_redeemr_trim4_binom.",name,".rds")


prepare_redeemr_v1<-function(redeemv_final, sample, reportpath){
    VariantsGTSummary<-redeemR.read(path=paste0(redeemv_final,sample),thr="S",Processed=T)
    redeemR<-Create_redeemR(VariantsGTSummary, maxctscut = 2)
    redeemR<- add_raw_fragment(redeemR)
    redeemR<-clean_redeem_removehot(redeemR)
    redeem_qc_report<-run_redeem_qc(redeemR,redeemR@HomoVariants)
    report_v1<-list(redeemR=redeemR,redeem_qc_report=redeem_qc_report)
 #  return(report_v1)
    saveRDS(report_v1,reportpath)
}

prepare_redeemr_trim4_v1 <- function(redeemv_final, sample, reportpath){
    VariantsGTSummary <- redeemR.read.trim(paste0(redeemv_final,sample), "S", Processed=T,
                                       "/VariantsGTSummary.S.trim4.RDS",edge_trim=4)
    redeemR<-Create_redeemR(VariantsGTSummary, maxctscut = 2)
    redeemR<- add_raw_fragment(redeemR)
    redeemR<-clean_redeem_removehot(redeemR)
    redeem_qc_report<-run_redeem_qc(redeemR,redeemR@HomoVariants)
    report_trim4_v1<-list(redeemR=redeemR,redeem_qc_report=redeem_qc_report)
    # return(report_trim4_v1)
    saveRDS(report_trim4_v1,reportpath)
    
}


prepare_redeemr_trim4_binom <- function(redeemv_final, sample, reportpath){
    VariantsGTSummary <- redeemR.read.trim(paste0(redeemv_final,sample), "S", Processed=T,
                                       "/VariantsGTSummary.S.trim4_binom.RDS",edge_trim=4)
    redeemR<-Create_redeemR_model(VariantsGTSummary)
    #redeemR<- clean_redeem(redeemR,fdr = 0.05)
    redeemR<-clean_redeem_removehot(redeemR)
    redeemR<- add_raw_fragment(redeemR,raw="RawGenotypes.Sensitive.StrandBalance")
    redeem_qc_report<-run_redeem_qc(redeemR, redeemR@HomoVariants)
    report_trim4_binom<-list(redeemR=redeemR,redeem_qc_report=redeem_qc_report)
#    return(report_trim4_binom)
    saveRDS(report_trim4_binom,reportpath)
}

## run
 prepare_redeemr_v1(redeemv_final,sample_folder,savename_1)
 prepare_redeemr_trim4_v1(redeemv_final,sample_folder,savename_2)
 prepare_redeemr_trim4_binom(redeemv_final,sample_folder,savename_3)

