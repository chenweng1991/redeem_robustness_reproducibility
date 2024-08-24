#library(ggtree)
#library(treeio)
#library(ape)
library("ggtree",lib.loc="/lab/solexa_weissman/cweng/Packages/miniconda/envs/r4.3_environment/lib/R/library")
library("treeio",lib.loc="/lab/solexa_weissman/cweng/Packages/miniconda/envs/r4.3_environment/lib/R/library")
library("ape",lib.loc="/lab/solexa_weissman/cweng/Packages/miniconda/envs/r4.3_environment/lib/R/library")
args <- commandArgs(trailingOnly = TRUE)
input=args[1]
output=args[2]

d<-readRDS(input)
tree<-as.treedata(nj(d))
current_directory <- getwd()
saveRDS(tree,paste(current_directory,output,sep="/"))

q