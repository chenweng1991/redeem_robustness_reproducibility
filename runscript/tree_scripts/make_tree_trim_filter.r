library(dplyr)
library(redeemR)
source("../../utils/utils.r")

dir1_redeem_binom <- "/lab/solexa_weissman/cweng/Projects/Collaborator/Caleb/github_redeem_plus/output/"
dir2="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/reproducibility_pub/data/redeemR_object_oldversion/"
dir3="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/reproducibility_pub/data/redeemV_final/"


args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
sample_folder <- args[2]

make_tree_trim <- function(redeem, raw1, raw2, trim=5, blacklist=c("Variants310TC", "Variants309CT", "Variants9979GA")){
    ## Step 1 make a trimed matrix list generate a list(Cts.Mtx,Cts.Mtx.bi) after trimming 
    cts.trim.mtx.bi<-make_ctxMatrix_edge0(redeem,raw1,raw2,"1","2",trim=trim)
    ## Step 2, remove variants with flag and black list
    cts.trim.remove.mtx.bi<-cts.trim.mtx.bi[,-which(colnames(cts.trim.mtx.bi) %in% c(blacklist))]
    cts.trim.remove.mtx.bi<-cts.trim.remove.mtx.bi[rowSums(cts.trim.remove.mtx.bi)>=2,]
    ## Step3, compute distance
    d.w_jaccard<-get_dist(cts.trim.remove.mtx.bi) %>% as.matrix
    ## Step4, make the tree
    tree<-as.treedata(nj(d.w_jaccard))
    list(tree=tree,cts.trim.remove.mtx.bi=cts.trim.remove.mtx.bi)
}


combined.redeem <- readRDS(paste0(dir2,"DN4_BMMC_HSPC_mitoTracing.sensitive_JACnj_Vassigned_p"))
sensitve.raw<-read.table(paste0(dir3,"Young1.T1.BMMC.Consensus.final","/RawGenotypes.Sensitive.StrandBalance"))
sensitve.raw.HPC <- read.table(paste0(dir3,"Young1.T1.HPC.Consensus.final","/RawGenotypes.Sensitive.StrandBalance"))  

result<- make_tree_trim(combined.redeem,sensitve.raw,sensitve.raw.HPC, trim=5)

saveRDS(result,"../../output/trees/Young1_BMMC_HSPC.trim5.tree.rds")