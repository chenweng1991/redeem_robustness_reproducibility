library(stringr)
library(Matrix)
library(Matrix.utils)
library(ggtree)
library(treeio)
library(ape)
#' run_redeem_qc
#'
#' This function generate qc plot to assess filtering strategy
#' @param redeem The redeemR object
#' @param homosets the homoplasmic mutation set
#' @param hotcall make sure hotcall mutations are not included
#' @export
run_redeem_qc <- function(redeem, homosets, hotcall= c("310_T_C","309_C_T","3109_T_C", "5764_C_T". "9979_G_A")){
    print("Make sure add_raw_fragment has been run  after clean_redeem")
    require(ggplot2)   
    require(gridExtra)
    require(igraph)    
    ## check the position distribution</span>
    p_pos<-subset(redeem@raw.fragment.uniqV, !variant %in% c(homosets,hotcall))  %>% ggplot()+aes(rel_position)+geom_histogram(bins=100,fill="white",color="black")+theme_cw1()+ggtitle("1+molecule")
    pos_1mol <- subset(redeem@raw.fragment.uniqV, !variant %in% c(homosets, hotcall)) %>% subset(., Freq==1)  %>% 
    ggplot()+aes(rel_position)+geom_histogram(bins=100,fill="white",color="black")+theme_cw1()+ggtitle("1 mol")
    grid.arrange(p_pos,pos_1mol)

    ## Add transition and transversion information
    redeem@V.fitered <- redeem@V.fitered %>% mutate(changes=add_changes(Variants)) %>% mutate(types=add_types(changes))
    redeem@raw.fragment.uniqV<-redeem@raw.fragment.uniqV %>% mutate(changes=add_changes(variant)) %>% mutate(types=add_types(changes)) 

    p_cell_maxcts <- ggplot(redeem@V.fitered) + aes(log2(CellN), log2(maxcts), color=types) +geom_point()+scale_color_brewer(palette="Set1")+theme_cw1()
    p_cell_meancts <- ggplot(redeem@V.fitered) + aes(log2(CellN), log2(PositiveMean_cts), color=types) +geom_point()+scale_color_brewer(palette="Set1")+theme_cw1()
    grid.arrange(p_cell_maxcts, p_cell_meancts, nrow=1)

    ## report transversion rate
    raw.fragment.uniqV_allhetero <- subset(redeem@raw.fragment.uniqV, !variant %in% homosets)
    raw.fragment.uniqV_1mol<-subset(redeem@raw.fragment.uniqV, Freq==1 & !variant %in% homosets)
    raw.fragment.uniqV_2mol<-subset(redeem@raw.fragment.uniqV, Freq==2 & !variant %in% homosets)
    raw.fragment.uniqV_3molplus<-subset(redeem@raw.fragment.uniqV, Freq>=3 & !variant %in% homosets)
    transversion_rate<-c(unweight=sum(redeem@V.fitered$types=="transversion")/nrow(redeem@V.fitered),
                         weight_freq_all=sum(raw.fragment.uniqV_allhetero$types=="transversion")/nrow(raw.fragment.uniqV_allhetero),
                         weight_freq_1mol=sum(raw.fragment.uniqV_1mol$types=="transversion")/nrow(raw.fragment.uniqV_1mol),
                         weight_freq_2mol=sum(raw.fragment.uniqV_2mol$types=="transversion")/nrow(raw.fragment.uniqV_2mol),
                         weight_freq_3molplus=sum(raw.fragment.uniqV_3molplus$types=="transversion")/nrow(raw.fragment.uniqV_3molplus))
    print(transversion_rate)

    GTsummary.filtered<- redeem@GTsummary.filtered %>% subset(.,!Variants %in% c(homosets,hotcall))

    ##check number of total mutations; mutation per cell, number of cells connected </span>

    Cts.Mtx <- dMcast(GTsummary.filtered,Cell~Variants,value.var = "Freq")
    Cts.Mtx.bi <- Cts.Mtx
    Cts.Mtx.bi[Cts.Mtx.bi>=1]<-1

    number_of_total_mut <- ncol(Cts.Mtx.bi)
    number_of_total_cells <- nrow(Cts.Mtx.bi)
    mut_per_cell <- rowSums(Cts.Mtx.bi)
    redeem_n0.adj.mtx<-CountOverlap_Adj(Cts.Mtx.bi,n=0)
    diag(redeem_n0.adj.mtx)<-0
    number_cell_connected<-rowSums(redeem_n0.adj.mtx) 
    mc <- graph_from_adjacency_matrix(redeem_n0.adj.mtx) %>% igraph::components(mode = "weak") 

    mc <- graph_from_adjacency_matrix(redeem_n0.adj.mtx) %>% igraph::components(mode = "weak") 
    fraction_in_component<-max(mc$csize)/nrow(redeem_n0.adj.mtx)

    print(paste0("number of cells: ", as.character(number_of_total_cells)))
    print(paste0("number of total mutations: ", as.character(number_of_total_mut)))
    print(paste0("number of mutations per cell: ", as.character(median(mut_per_cell))))
    print(paste0("number of cells connected: ", as.character(median(number_cell_connected))))
    print(paste0("fraction of cells in component: ", as.character(median(fraction_in_component))))

    report_metric<-list(number_of_total_mut=number_of_total_mut, 
                        number_of_total_cells=number_of_total_cells,
                        mut_per_cell=mut_per_cell, 
                        fraction_in_component=fraction_in_component,
                        number_cell_connected=number_cell_connected )

    plots <- list(p_pos=p_pos,pos_1mol=pos_1mol, p_cell_maxcts=p_cell_maxcts,p_cell_meancts=p_cell_meancts)

    return(list(plots=plots, transversion_rate=transversion_rate, report_metric=report_metric))
}




## Functions added (included in redeemR)  2024-8-5

## Add frequency to raw fragments
add_freq_raw<-function(raw){
GiveName <- c("UMI", "Cell", "Pos", "Variants", "Call", "Ref", "FamSize", "GT_Cts", "CSS", "DB_Cts", "SG_Cts", "Plus", "Minus", "Depth")
colnames(raw)<-GiveName
raw$CellVar<-paste(raw$Cell,raw$Variants,sep="_")
raw.gtsummary<-GTSummary(raw)
raw<-merge(raw,raw.gtsummary[,c("Var1","Freq")],by.x="CellVar",by.y="Var1",all.x = T)
return(raw)
}

#' Function needed to compute the 
make_position_df_3.4<-function(in_df){
    first <- str_split_fixed(in_df[,"UMI"], "_", 3)[, c(2)] %>% 
        as.numeric()
    last <- str_split_fixed(in_df[,"UMI"], "_", 3)[, c(3)] %>% 
        as.numeric()
    start <- pmin(first, last)
    end <- pmax(first, last)
    df <- data.frame(UMI=in_df[,"UMI"],start = start, end = end, pos = in_df$Pos, 
        variant = in_df$Variants) %>% mutate(length = end - start) %>% 
        mutate(rel_position = (pos - start)/length,Freq=in_df[,"Freq"])
    df$edge_dist<-pmin(abs(df$pos-df$start),abs(df$pos-df$end))
    return(df)
}

#' make_position, but remain all the raw genotyp information
#' Input is the ind_df or 
#' 
#' @param in_df  raw genotype by redeemV
#' @export
## 
make_position_df_3.5<-function(in_df){
    first <- str_split_fixed(in_df[,"UMI"], "_", 3)[, c(2)] %>% 
        as.numeric()
    last <- str_split_fixed(in_df[,"UMI"], "_", 3)[, c(3)] %>% 
        as.numeric()
    start <- pmin(first, last)
    end <- pmax(first, last)
    df <- data.frame(in_df,start = start, end = end, pos = in_df$Pos, 
        variant = in_df$Variants) %>% mutate(length = end - start) %>% 
        mutate(rel_position = (pos - start)/length,Freq=in_df[,"Freq"])
    df$edge_dist<-pmin(abs(df$pos-df$start),abs(df$pos-df$end))
    return(df)
}



#' Produce a raw fragment table with frequency (how many cells) and the reletive distance
#' from redeemR object, 
#' 
#' @param redeemR  a redeemR object 
#' @export
add_raw_fragment <- function(redeemR,raw="RawGenotypes.Sensitive.StrandBalance"){
    if ("edge_trim" %in% names(redeemR@para)){
        edge_trim <- as.numeric(redeemR@para["edge_trim"])
    }else{
        edge_trim <- 0
    }
    print(paste0("It has benn edge trimmed by ", as.character(edge_trim), " bp"))
    redeemR.raw <- read.table(paste(redeemR@attr$path,raw,sep="/"))
    filtered.variants <- unique(redeemR@GTsummary.filtered$Variants)
    redeemR.raw.passfilter <- subset(redeemR.raw, V4 %in% filtered.variants)
    raw.pos<- redeemR.raw.passfilter %>% add_freq_raw() %>% make_position_df_3.4() %>% filter(edge_dist>=edge_trim)
    redeemR@raw.fragment.uniqV <- raw.pos
    return(redeemR)
}


#' Convinient function that takes raw fragment, and out put fragment with frequency
#' Produce KS test results, input is the redeem object
#' from redeemR object, 
#' 
#' @param raw
#' @export
add_freq_raw<-function(raw){
GiveName <- c("UMI", "Cell", "Pos", "Variants", "Call", "Ref", "FamSize", "GT_Cts", "CSS", "DB_Cts", "SG_Cts", "Plus", "Minus", "Depth")
colnames(raw)<-GiveName
raw$CellVar<-paste(raw$Cell,raw$Variants,sep="_")
raw.gtsummary<-GTSummary(raw)
raw<-merge(raw,raw.gtsummary[,c("Var1","Freq")],by.x="CellVar",by.y="Var1",all.x = T)
return(raw)
}



#' Produce KS test results, input is the redeem object
#' from redeemR object, 
#' 
#' @param redeemR  a redeemR object 
#' @export
make_ks_test_df <- function(redeemR){
print("Make sure redeemR@raw.fragment.uniqV exist, produced by redeemR@raw.fragment.uniqV <- add_raw_fragment(redeemR)")
  redeemR@raw.fragment.uniqV %>%
  nest(data = c(-variant)) %>%
    mutate(D_stat=purrr::map(data, ~(ks.test(.x$rel_position, "punif",0,1))),
           tidied = purrr::map(D_stat, broom::tidy)) %>%
    mutate(n = map_dbl(data, nrow)) %>%
    tidyr::unnest(tidied) %>%
    dplyr::select(variant, statistic, p.value, n)
}

######
###### Functions added 24-08-09

## Internal function for Create_redeemR_model
goodness_of_fit_test <- function(trials,coverage = 30) {
  # Number of trials
  n_trials <- length(trials)
  # Handle possible coverage exceed
  if (max(trials)>coverage){
    coverage <- max(trials)
  }  
  # Observed frequencies
  observed_frequencies <- table(factor(trials, levels = 0:coverage))
  # Single coin model
  total_heads <- sum(trials)
  total_flips <- n_trials * coverage
  p_hat <- total_heads / total_flips
  # Expected frequencies under the single binomial model
  expected_frequencies <- n_trials * dbinom(0:coverage, size = coverage, prob = p_hat)
  # Handle the 0 in expected frequency
  expected_frequencies[expected_frequencies==0]<- 1e-318  
  # Chi-squared test statistic
  chi2_stat <- sum((observed_frequencies - expected_frequencies)^2 / expected_frequencies)
  # p-value from chi-squared distribution with coverage degrees of freedom (coverage + 1 possible outcomes - 1 estimated parameter)
  p_value <- pchisq(chi2_stat, df = coverage, lower.tail = FALSE)
  return(list(chi2_stat = chi2_stat, p_value = p_value))
}

## Internal function, executed by Create_redeemR_model
run_binomial_noise_removal <- function(redeem){
    require(qvalue)
    if ("chi" %in% colnames(redeem@V.fitered)){
        stop("goodness_of_fit_test has been run on this dataset")
    }else{    
        Mtx <- redeem@Cts.Mtx
        stats <- c()
        pvalues <- c()
        for (i in 1:ncol(Mtx)){
            v<- as.integer(redeem@Cts.Mtx[,i])
            pos <- as.numeric(gsub("Variants([0-9]+)[A-Za-z]{2}$", "\\1", colnames(Mtx)))
            cov<-as.integer(as.numeric(redeem@DepthSummary$Pos.MeanCov[pos[i],"meanCov"]))
            res<-goodness_of_fit_test(v,cov)
            stats<-c(stats,res$chi2_stat)
            pvalues<-c(pvalues,res$p_value)
        }
        qvalues <- qvalue(pvalues)$qvalues
        redeem@V.fitered<-merge(redeem@V.fitered,data.frame(Variants=convert_variant(colnames(Mtx)),pvalues=pvalues, qvalues=qvalues, chi=stats),all = T)
        return(redeem)
        }    
}

#' Create_redeemR_model
#'
#' This function is to create redeemR with basic information, modeled by binomial distribution
#' maxctscut is removed, replaced by the modeling, will generate statistics and p value
#' @param VariantsGTSummary simply put GTSummary (Generated by redeemR.read) 
#' @param qualifiedCellCut The minimum median mitochondrial coverage for a qualified cell, default is 10
#' @param OnlyHetero If only consider the heteroplasmy variants, default is T
#' @param VAFcut only use variants with VAF smaller than VAFcut. Default is 1.  We can use smaller value to constrain into only using rare variants
#' @param Cellcut only use variants with at least cellcut cells carry
#' @return redeemR class
#' @export
#' @import Seurat ape phytools phangorn tidytree ggtreeExtra
#' @importFrom ggtree ggtree
Create_redeemR_model<-function(VariantsGTSummary=VariantsGTSummary,qualifiedCellCut=10,VAFcut=1,Cellcut=2){
 if ("edge_trim" %in% names(attributes(VariantsGTSummary))){
        edge_trim <- as.numeric(attr(VariantsGTSummary,"edge_trim"))
    }else{
        edge_trim <- 0
    }
CellMeta<-subset(attr(VariantsGTSummary,"depth")[["Cell.MeanCov"]],meanCov>=qualifiedCellCut)
names(CellMeta)[1]<-"Cell"
VariantsGTSummary.feature<-Vfilter_v4(VariantsGTSummary,Min_Cells = Cellcut, Max_Count_perCell = 1, QualifyCellCut = qualifiedCellCut)
GTsummary.filtered<-subset(VariantsGTSummary,Variants %in% VariantsGTSummary.feature$Variants & Cell %in% CellMeta$Cell)
ob<-new("redeemR")
ob@GTsummary.filtered<-GTsummary.filtered
ob@CellMeta<-CellMeta
ob@V.fitered=VariantsGTSummary.feature
ob@HomoVariants<-attr(VariantsGTSummary.feature,"HomoVariants")
ob@UniqueV<-VariantsGTSummary.feature$Variants
ob@DepthSummary<-attr(VariantsGTSummary,"depth")
ob@para<-c(Threhold=attr(VariantsGTSummary,"thr"),qualifiedCellCut=qualifiedCellCut,VAFcut=VAFcut,Cellcut=Cellcut,edge_trim=edge_trim)
ob@attr<-list(path=attr(VariantsGTSummary,"path"))
ob<-Make_matrix(ob,onlyhetero=T)
ob<-run_binomial_noise_removal(ob)
return(ob)
}



#' clean_redeem
#'
#' This function is to clean redeem by filtering both V.fitered and GTsummary.filtered by qvalues
#' 
#' @param ob redeem object
#' @param fdr fdr cutoff, default is 0.05
#' @export
clean_redeem <-function(ob,fdr=0.05){
    ob@V.fitered <- subset(ob@V.fitered,qvalues<=0.05)
    ob@GTsummary.filtered<-subset(ob@GTsummary.filtered, Variants %in% ob@V.fitered$Variants)
    ob<-Make_matrix(ob,onlyhetero=T)
    ob@UniqueV <- ob@V.fitered$Variants
    return(ob)
}


#' clean_redeem_removehotcall
#'
#' This function is to clean redeem by filtering both V.fitered and GTsummary.filtered by qvalues
#' 
#' @param ob redeem object
#' @param hotcall fdr cutoff, default is 0.05
#' @export
clean_redeem_removehot <-function(ob,hotcall= c("310_T_C","9979_G_A","3109_T_C")){
    ob@V.fitered <- subset(ob@V.fitered,!Variants %in% hotcall)
    ob@GTsummary.filtered<-subset(ob@GTsummary.filtered, !Variants %in% hotcall)
    ob<-Make_matrix(ob,onlyhetero=T)
    ob@UniqueV <- ob@V.fitered$Variants
    return(ob)
}

######
###### general helper Functions added 24-08-09

## Internal function to convert the variant names, implemented by convert_variant
convert_variant_1 <- function(input_string) {
  if (grepl("^Variants[0-9]+[A-Za-z]{2}$", input_string)) {
    # Convert from "Variants10000GA" to "10000_G_A"
    number <- gsub("Variants([0-9]+)[A-Za-z]{2}$", "\\1", input_string)
    letters <- gsub("^Variants[0-9]+([A-Za-z]{2})$", "\\1", input_string)
    output_string <- paste0(number, "_", substr(letters, 1, 1), "_", substr(letters, 2, 2))
  } else if (grepl("^[0-9]+_[A-Za-z]_[A-Za-z]$", input_string)) {
    # Convert from "10000_G_A" to "Variants10000GA"
    parts <- strsplit(input_string, "_")[[1]]
    output_string <- paste0("Variants", parts[1], parts[2], parts[3])
  } else {
    stop("Input string format is not recognized.")
  }
  return(output_string)
}

#' convert_variant
#'
#' @param x this is a vector of variant names, either way for example from '10000_G_A' to/from 'Variants10000GA'
#' @return a vector of strings 
#' @export
convert_variant <- function(x){
    res<-sapply(x,convert_variant_1)
    return(as.character(res))
}

#' add_changes
#'
#' @param variant given a variant, output the changes 
#' @return changes
#' @export
add_changes <-function(variant){
    changes<-sub("^\\d+_", "",variant)
    return(changes)
}

#' add_types
#'
#' @param changes given a changes, output the type (transition or transversion)
#' @return changes
#' @export
add_types <- function(changes){
    types<- ifelse(changes %in% c("C_T","G_A","T_C","A_G"),"transition","transversion")
    return(types)
}

#' Annotate_base_change
#' 
#' @param  input redeem object,  it takes V.fitered, add the nucleotide change, 
#' @return redeem object with the V.fitered slot modified
Annotate_base_change <- function(redeem){
    redeem@V.fitered<- redeem@V.fitered %>% mutate(changes=add_changes(Variants)) %>% mutate(types=add_types(changes))
    return(redeem)
}


#' Define a custom theme function
#' 
#' @param  axis_title_size
#' @param  axis_text_size
#' @export
theme_cw1 <- function(axis_title_size = 20, axis_text_size = 15) {
    theme_classic()+theme(
    axis.title = element_text(size = axis_title_size,color="black"),  # Set axis title size
    axis.text = element_text(size = axis_text_size,color="black")     # Set axis text size
  )
}

#' Define a custom theme function
#' 
#' @param  axis_title_size
#' @param  axis_text_size
#' @export
theme_cw2 <- function(axis_title_size = 20, axis_text_size = 15) {
    theme_classic()+theme(
    axis.title = element_text(size = axis_title_size,color="black"),  # Set axis title size
    axis.text = element_text(size = axis_text_size,color="black")     # Set axis text size
  )
}


#' CountOverlap_Adj
#' function to count the connectedness (adjacency matrix), or the number of cells sharing more than n variants with the given cell
#' @param  M
#' @param  n 
#' @export
CountOverlap_Adj<-function(M,n=0){
    require(Matrix)
    # Total <- Matrix::rowSums(M)
    a <- M %*% Matrix::t(M)
    a[a<=n]<-0
    a[a>n]<-1
    return(a)
}


#' call_lmhc
#' Legacy function to call LMHC, in the latest redeem, the mean count is included in the the V.fitered, and LMHC can be directly filtered by 
#' @param  M
#' @param  n 
#' @export
call_lmhc<-function(redeemR,cuoff_meancount=1.3,cutoff_celln=20){
meancounts<-make_meancount(redeemR)
LMHC<-subset(meancounts,meancounts<cuoff_meancount & ncells>cutoff_celln) %>% row.names
return(LMHC)
}

#' make_meancount
#' Internal function to compute the meancount
#' This is a legacy function, in the latest redeem, the mean count is included in the the V.fitered
#' @param  redeemR object
make_meancount<-function(redeemR){
meancounts<-redeemR@Cts.Mtx %>% as.matrix %>% apply(.,2,function(x){mean(x[x>0])}) %>% data.frame(meancounts=.)
ncells<-redeemR@Cts.Mtx %>% as.matrix %>% apply(.,2,function(x){sum(x>0)}) %>% data.frame(ncells=.)
meancounts$variants<-gsub("Variants(\\d+)([A-Z])([A-Z])", "\\1_\\2_\\3", row.names(meancounts))
meancounts<-Tomerge_v2(meancounts,ncells)
return(meancounts)
}

get_dist<-function(Cts.Mtx.bi){
    V.weight<-data.frame(weight=1-CellPCT$muRate)
    V.weight$Variants<-paste("Variants",gsub("_","",CellPCT$Variant),sep="")
    weight<-data.frame(Variants=colnames(Cts.Mtx.bi)) %>% merge(.,V.weight,by="Variants",all.x = T,sort = F) %>% .$weight
    weight[is.na(weight)]<-1
    d.w_jaccard<-quick_w_jaccard(Cts.Mtx.bi,w=weight)
    return(d.w_jaccard)
}

get_hyper<- function(cut=0.01){
    data(CellPCT)
    Hypermutatble<-subset(CellPCT,(PCT.D4BM >cut | PCT.D4HPC>cut | PCT.D4HSC >cut) &  (PCT.D9BM>cut |PCT.D9HPC>cut | PCT.D9HSC>cut) & (PCT.D1BM>cut | PCT.D1HPC>cut ) )$Variants
    return(Hypermutatble)
}

get_hype_v2<- function(cut=0.01){
    Hypermutatble<- subset(CellPCT.update, (PCT.D4BM >cut | PCT.D4HPC>cut | PCT.D4HSC >cut) & (PCT.D9BM > cut | PCT.D9HPC > 
        cut | PCT.D9HSC > cut) & (PCT.D1BM > cut | PCT.D1HPC > 
        cut) & (Old1_BM > cut | Old1_HSPC> cut) & (Old2_BM > cut | Old2_HSPC> cut))$Variants
    return(Hypermutatble)
}




#' Function to read in multiple runs of redeemV outputs with edge trimming, default is trimming 5bp
#'
#' This function allows you to read raw data from multiple XX/final folder, the output from redeemV
#' It process the data same way as CW_mgatk.read but need to specify one threadhold(thr)
#' @param paths a vector of multiple path strings The XX/final folder, the output from mitoV
#' @param names The names for the multiple input, same order with folder
#' @param suffix The suffix for the multiple input, which is added at the end of cell name, separated by "."
#' @param edge_trim  how many bp to be trimmed, default is 5, 
#' @return VariantsGTSummary  combined VariantsGTSummary
#' @export
redeemR.read.multiple.trim<-function(paths,names, suffix, edge_trim=5){
    print("This version hardcode in using sensitive or S, if other consensus is needed, feel free to modify this source code by changing RawGenotypes.Sensitive.StrandBalance")
    thr="S"
    GiveName<-c("UMI","Cell","Pos","Variants","Call","Ref","FamSize","GT_Cts","CSS","DB_Cts","SG_Cts","Plus","Minus","Depth")
    ## Add suffix to cell names and concat the RawGenotypes
    RawGenotypes.concat <- c()
    for (i in 1: length(paths)){
        RawGenotypes<-read.table(paste(paths[i],"/RawGenotypes.Sensitive.StrandBalance",sep=""))
        RawGenotypes$V2 <- paste(RawGenotypes$V2, suffix[i], sep=".")
        RawGenotypes.concat <- rbind(RawGenotypes.concat,RawGenotypes)
    }    
    ## Add the column name
    colnames(RawGenotypes.concat)<-GiveName
    ## Annotate the rawGenotype by dstance, etc
    RawGenotypes.annotated <- RawGenotypes.concat %>% add_freq_raw() %>% make_position_df_3.5()
    ## Trim
    RawGenotypes.trimed<- filter(RawGenotypes.annotated,edge_dist>=edge_trim)
    ## Do GTSummary, this part is the same as regular
    VariantsGTSummary.trimed<-GTSummary(RawGenotypes.trimed)
    VariantsGTSummary<-GTSummary(RawGenotypes.concat)
    VariantsGTSummary<- VariantsGTSummary[,c("Var1","Freq")] %>% rename (Freq_before_trim=Freq) %>% 
    merge(VariantsGTSummary.trimed,.,by="Var1") %>% mutate(depth=depth-(Freq_before_trim-Freq))
    VariantsGTSummary$hetero<-with(VariantsGTSummary,Freq/depth)
    ## Compute for the depth, this part integrate and compute the combined depth
    Pos.MeanCov<-data.frame(pos=1:16569)
    Cell.MeanCov<-c()
    cell.numbers<-c()
    for (i in 1:length(paths)){
        depth_i<-DepthSummary(paths[i])
        pos.info <- depth_i[[1]]
        cell.info <- depth_i[[2]]
        cell.info$V1 <- paste(cell.info$V1, suffix[i], sep=".")
        Pos.MeanCov<-cbind(Pos.MeanCov,pos.info[,"meanCov",drop=F]) 
        Cell.MeanCov<-rbind(Cell.MeanCov,cell.info)
        cell.numbers<-c(cell.numbers,nrow(cell.info))
    }
    Pos.MeanCov.final<- data.frame(V2=1:16569, meanCov=rowSums(Pos.MeanCov[,2:ncol(Pos.MeanCov)] * cell.numbers)/sum(cell.numbers))  ## Compute final cov by pos
    ## add attributes
    attr(VariantsGTSummary,"depth")<-list(Pos.MeanCov=Pos.MeanCov.final, Cell.MeanCov=Cell.MeanCov)
    attr(VariantsGTSummary,"thr")<-thr
    attr(VariantsGTSummary,"path")<-paths
    attr(VariantsGTSummary,"edge_trim")<-edge_trim
    attr(VariantsGTSummary,"combined")<-names
    attr(VariantsGTSummary,"suffix")<-suffix
    return(VariantsGTSummary)
}



#' add_mutation_type
#' add a  
#' 
#' @param redeemR  a redeemR object 
#' @export
add_mutation_type <- function(redeem){
   redeem@V.fitered$types<- redeem@V.fitered$Variants %>% add_changes() %>% add_types()
   return(redeem) 
}


#### 2024-8-11, added for Notebook-01_edge_quantification
