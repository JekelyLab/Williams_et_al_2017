

'''
Author: Olivier Mirabeau
Date: 03/12/2017
This script is designed to take two Single cell RNA count matrices (here single cells datasets from the Romanov et al. 2017 PMID: 27991900 and
Campbell et al. 2017 PMID: 28166221 studies and the Platynereis dataset from Williams et al. 2017)
from two different species with an orthology structure described in the variable pdum.TF.genes.ortho.list
and derive correlation plots of all vs all markers, with a conserved order between the orthologous markers sets
'''

###################################################################
############## setting root directory #############################
###################################################################

setwd("E:/Williams_Elife2017")

###################################################################
############## loading of libraries ############################### 
###################################################################

library(minerva)
library(corrplot)
library(stringr)

##########################################################################
############## loading of the Romanov et al. (2017) dataset ############## 
##########################################################################

load(verbose=T, "mouse_datasets/romanov2017/GSE74672_expressed_molecules.RData")
mouse.names.romanov <-  as.vector(raw.romanov.hypo.cells[,1])
mouse.romanov.dat <- do.call(cbind, lapply(1:ncol(raw.romanov.hypo.cells[,-1]), function(x) as.numeric(as.character(raw.romanov.hypo.cells[,-1][,x]))))
mouse.romanov.dat <- log10(mouse.romanov.dat+1)
rownames(mouse.romanov.dat) <- mouse.names.romanov
colnames(mouse.romanov.dat) <- colnames(raw.romanov.hypo.cells[,-1])

###########################################################################
############## loading of the Campbell et al. (2017) dataset ############## 
###########################################################################

load(verbose=T, "mouse_datasets/campbell2017/GSE93374_Merged_all_020816_DGE.RData")

campbell.2017.metadata <- read.table(header=T, sep="\t", "mouse_datasets/campbell2017/GSE93374_cell_metadata.txt")
campbell.2017.sc.table.neurons <- campbell.2017.sc.table[ , which(as.vector(campbell.2017.metadata$X8.clust_all_neurons)=="neuron")]
groups.X2 <- as.vector(unlist(campbell.2017.metadata["X2.group"]))[which(as.vector(campbell.2017.metadata$X8.clust_all_neurons)=="neuron")]
campbell.2017.sc.neurons.dat <- do.call(cbind, lapply(1:ncol(campbell.2017.sc.table.neurons), function(x) as.numeric(as.character(campbell.2017.sc.table.neurons[,x]))))
rownames(campbell.2017.sc.neurons.dat) <- rownames(campbell.2017.sc.table.neurons)
colnames(campbell.2017.sc.neurons.dat) <- colnames(campbell.2017.sc.table.neurons)
riken.ids.campbell <- rownames(campbell.2017.sc.neurons.dat)

###########################################################################
############## loading of the Williams et al. (2017) dataset ############## 
###########################################################################

pdum.sc.table <- read.table(sep="\t", "Platynereis_dataset/single_cells_RNA-SEQ_RPKM-normalize-log10-clustered_joined_LR-cells.csv")
pdum.sc.dat <- do.call(cbind, lapply(1:ncol(pdum.sc.table[-1,-c(1:3)]), function(x) as.numeric(as.character(pdum.sc.table[-1,-c(1:3)][,x]))))
rownames(pdum.sc.dat) <- as.vector(pdum.sc.table[-1,1])
colnames(pdum.sc.dat) <- as.vector(unlist(pdum.sc.table[1,-c(1:3) ]))

#############################################  
############### parameters ##################
#############################################

pdum.or.mouse <- "mouse"  # c("pdum", "mouse")
type.of.mouse.data <- "romanov" # c("romanov", "campbell")
limit.pval.corrected <- 0.05
randomisation_folder.mouse <- "randomisation_romanov" ## for the romanov et al dataset

############### parameters #############

# from now on we work either on the Romanov dataset

annots.mouse <- rownames(mouse.romanov.dat)

####################################################################
############## loading the orthology table #########################
####################################################################

ortho.table.pdum.mouse <- read.table("ortholog_table_pdum_mouse_clean.txt", sep="\t", header=F)
mouse.TF.genes.ortho.list <- lapply(1:nrow(ortho.table.pdum.mouse),function(x) ortho.table.pdum.mouse[x,-1][ortho.table.pdum.mouse[x,-1]!=""])
match.TF.names.mouse <- lapply(1:length(mouse.TF.genes.ortho.list), function(x) match(mouse.TF.genes.ortho.list[[x]], annots.mouse))
pdum.TF.genes.ortho.list <- lapply(1:nrow(ortho.table.pdum.mouse),function(x) as.vector(ortho.table.pdum.mouse[x,1][ortho.table.pdum.mouse[x,1]!=""]))



############### parameters #############
pdum.or.mouse<-"mouse"  # c("pdum", "mouse")
pval.correction <- TRUE
limit.pval.uncorrected <- 1e-3
limit.pval.corrected <- 0.05

########################################################################################      
##### we will have to run this at least twice, once for the platynereis datset, once for the mouse datset  #########################
########################################################################################    

annots <- NULL
tab.current <- NULL
if (pdum.or.mouse=="mouse"){
    annots <- mouse.names.romanov
    tab.current <- mouse.romanov.dat
} else if (pdum.or.mouse=="pdum"){
  annots <- rownames(pdum.sc.dat)
  tab.current <- pdum.sc.dat
}
names.sym <- annots

if (is.null(colnames(tab.current))){
  colnames(tab.current) <- as.character(1:ncol(tab.current))
}


regex.frame <- "frame[-]{0,1}\\d_ORF_\\d+_(.*)"
needs.names.extraction <- grepl("frame", pdum.TF.genes.ortho.list)
pdum.TF.genes.ortho.list.extracted <- unlist(pdum.TF.genes.ortho.list)
pdum.TF.genes.ortho.list.extracted[which(needs.names.extraction)] <- str_match(pdum.TF.genes.ortho.list.extracted, regex.frame)[,2][which(needs.names.extraction)]
manual.annots.ortho.pdum.table <- read.table("Platynereis_dataset/Pdum_sequences_liz_orthology_namesfixed.txt", sep="\t", header=T)
does.not.need.contignames <- grepl("Pdu_", pdum.TF.genes.ortho.list.extracted)  # doesn't apply for the peptide genes!
pdum.TF.contignames.ortho.list.extracted <- trimws(as.vector(manual.annots.ortho.pdum.table[match(pdum.TF.genes.ortho.list.extracted, as.vector(manual.annots.ortho.pdum.table[,1])), 2]))
pdum.TF.contignames.ortho.list.extracted[does.not.need.contignames] <- pdum.TF.genes.ortho.list.extracted[does.not.need.contignames]

##################################################     
##### we do all vs. all  #########################
##################################################     

list.1 <- NULL
list.2 <- NULL
if (pdum.or.mouse=="mouse"){
  list.1 <- mouse.TF.genes.ortho.list  
  list.2 <- mouse.TF.genes.ortho.list 
} else if (pdum.or.mouse=="pdum"){
  list.1 <- lapply(1:length(pdum.TF.contignames.ortho.list.extracted), function(x) pdum.TF.contignames.ortho.list.extracted[x])
  list.2 <- lapply(1:length(pdum.TF.contignames.ortho.list.extracted), function(x) pdum.TF.contignames.ortho.list.extracted[x])
}

############################################################
##### initialisation of variables ##########################
############################################################

peprec.corr.mat <- matrix(0,length(list.1), length(list.2))
peprec.pval.mat <- matrix(0,length(list.1), length(list.2))
peprec.or.mat <- matrix(0,length(list.1), length(list.2))
peprec.pval.fisher.mat <- matrix(0,length(list.1), length(list.2))
peprec.pep.mat <- matrix(0,length(list.1), length(list.2))
peprec.rec.mat <- matrix(0,length(list.1), length(list.2))
peprec.numintersect.mat <- matrix(0,length(list.1), length(list.2)) # number of common cells expressing more than 1 read of both genes
peprec.numcellnon0.1.mat <- matrix(0,length(list.1), length(list.2))
peprec.numcellnon0.2.mat <- matrix(0,length(list.1), length(list.2))
peprec.mic.mat <- matrix(0,length(list.1), length(list.2))
peprec.besti.mat <- matrix(0,length(list.1), length(list.2)) 
peprec.bestj.mat <- matrix(0,length(list.1), length(list.2)) 

############################################################
########### calculation of correlations - main loop ########
############################################################

for (x in 1:length(list.1)) {
  pep.vec <- unlist(unname(list.1[[x]]))
  bestpep <- pep.vec[[1]]
  for (y in 1:length(list.2)) {
    rec.vec <- unlist(unname(list.2[[y]]))
    print(paste0("num of PepSystems: ",x, "  ", y))
    min.pval <- 2 
    max.coexpr <- 0 
    bestrec <- rec.vec[[1]]
    min.pval.ctest <- 1
    if (type.of.test=="fisher"){
      best.numinter <- 0
      best.numcellnon0.1 <- 0
      best.numcellnon0.2 <- 0
      min.or.fisher <- -1000
      min.pval.fisher <- 1
    }
    bestrec <- rec.vec[[1]]
    besti <- 1
    bestj <- 1
    for (i in 1:length(pep.vec)){
      for (j in 1:length(rec.vec)){
        ind.pep <- match(pep.vec[i], names.sym)
        ind.rec <- match(rec.vec[j], names.sym)
        if (length(ind.pep) & length(ind.rec)>0){
          if (!is.na(ind.pep) & !is.na(ind.rec)){
            if (!(all(tab.current[ind.pep, ]==0) | all( tab.current[ind.rec, ]==0))){
              coexpr <- cor(tab.current[ind.pep, ], tab.current[ind.rec, ])
              # test of correlation
              ctest <- cor.test(tab.current[ind.pep, ], tab.current[ind.rec, ])
              coexpr <- ctest[[4]]
              pval.ctest <- ctest[[3]]
              # Fisher test on non-zero cells
              if (type.of.test=="fisher"){
                non0.pep <- names(which(tab.current[ind.pep, ]!=0))
                non0.rec <- names(which(tab.current[ind.rec, ]!=0))
                all.cell.names <- colnames(tab.current)
                ftest <- fisher.test.2genelists(non0.pep, non0.rec, all.cell.names)
                num.intersect <- length(intersect(non0.pep, non0.rec))
                ftest <- fisher.test.numeric(length(non0.pep), length(non0.rec), num.intersect, length(all.cell.names))
                or.fisher <- ftest[["estimate"]]
                pval.fisher <- ftest[["p.value"]]
              }
              
              pval <- pval.ctest
              test.used <- (pval < min.pval)
              pep <- pep.vec[i]
              rec <- rec.vec[j]
              if (test.used){
                max.coexpr <- coexpr
                min.pval <- pval
                min.pval.ctest <- pval.ctest
                bestpep <- pep
                bestrec <- rec
                if (type.of.test=="fisher"){
                  best.numinter <- num.intersect
                  best.numcellnon0.1 <- length(non0.pep)
                  best.numcellnon0.2 <- length(non0.rec)
                  min.pval.fisher <- pval.fisher
                  min.or.fisher <- or.fisher
                }
                besti <- i
                bestj <- j
                
              }
            }
            
          } 
        } 
      }
    }
    peprec.corr.mat[x,y] <- coexpr
    peprec.pval.mat[x,y] <- pval.ctest
    
    peprec.pep.mat[x,y] <- bestpep
    peprec.rec.mat[x,y] <- bestrec
    if (type.of.test=="fisher"){
      peprec.numintersect.mat[x,y] <- best.numinter
      peprec.numcellnon0.1.mat[x,y] <- best.numcellnon0.1
      peprec.numcellnon0.2.mat[x,y] <- best.numcellnon0.2
      peprec.pval.fisher.mat[x,y] <- pval.fisher
      peprec.or.mat[x,y] <- or.fisher
    }
    #   peprec.mic.mat[x,y] <- bestmic
    peprec.besti.mat[x,y] <- besti
    peprec.bestj.mat[x,y] <- bestj
  }
}


###############################################################
########### calculation of MICs on only the best pairs ########
###############################################################

for (x in 1:length(list.1)) {
  pep.vec <- unlist(unname(list.1[[x]]))
  bestpep <- pep.vec[[1]]
  for (y in 1:length(list.2)) {
    rec.vec <- unlist(unname(list.2[[y]]))
    ind.pep <- match(pep.vec[peprec.besti.mat[x,y]], names.sym)
    ind.rec <- match(rec.vec[peprec.bestj.mat[x,y]], names.sym)
    mic <- 0
    num.intersect <- 0
    non0.pep <- 0
    non0.rec <- 0
    or.fisher <- -1000
    pval.fisher <- 1
    if (length(ind.pep) & length(ind.rec)>0){
      if (!is.na(ind.pep) & !is.na(ind.rec)){
        if (!(all(tab.current[ind.pep, ]==0) | all( tab.current[ind.rec, ]==0))){ 
          coexpr <- cor(tab.current[ind.pep, ], tab.current[ind.rec, ])
          mic.all <- mine(tab.current[ind.pep, ],tab.current[ind.rec, ], master=TRUE, use="all.obs")
          mic = mic.all$MIC
          
          # Fisher test on non-zero cells
          non0.pep <- names(which(tab.current[ind.pep, ]!=0))
          non0.rec <- names(which(tab.current[ind.rec, ]!=0))
          all.cell.names <- colnames(tab.current)
          num.intersect <- length(intersect(non0.pep, non0.rec))
          ftest <- fisher.test.numeric(length(non0.pep), length(non0.rec), num.intersect, length(all.cell.names))
          or.fisher <- ftest[["estimate"]]
          pval.fisher <- ftest[["p.value"]]
          
        }
      }
    }
    peprec.mic.mat[x,y] <- mic
    peprec.or.mat[x,y] <- or.fisher
    peprec.pval.fisher.mat[x,y] <- pval.fisher
    peprec.numintersect.mat[x,y] <- num.intersect
    peprec.numcellnon0.1.mat[x,y] <- length(non0.pep)
    peprec.numcellnon0.2.mat[x,y] <- length(non0.rec)
  }
}


###########################################
########### naming considerations #########
###########################################

if (pdum.or.mouse=="mouse"){
  pep.backtocommonname <- FALSE  
  rec.backtocommonname <- FALSE 
} else if (pdum.or.mouse=="pdum"){
  pep.backtocommonname <- TRUE ## 
  rec.backtocommonname <- TRUE ## 
}

peprec.pep.mat.commonname <- peprec.pep.mat
if (pep.backtocommonname) {
  for (i in 1:nrow(peprec.pep.mat)) {
    for (j in 1:ncol(peprec.pep.mat)) {
      if (!grepl("Pdu_", peprec.pep.mat[i,j])){
        peprec.pep.mat.commonname[i,j] <- as.vector(manual.annots.ortho.pdum.table[match(peprec.pep.mat[i,j], as.vector(manual.annots.ortho.pdum.table[,2])), 1])
      }
    }
  }
}
peprec.rec.mat.commonname <- peprec.rec.mat
if (rec.backtocommonname) {
  for (i in 1:nrow(peprec.rec.mat)) {
    for (j in 1:ncol(peprec.rec.mat)) {
      if (!grepl("Pdu_", peprec.rec.mat[i,j])){
        peprec.rec.mat.commonname[i,j] <- as.vector(manual.annots.ortho.pdum.table[match(peprec.rec.mat[i,j], as.vector(manual.annots.ortho.pdum.table[,2])), 1])
      }
    }
  }
}

########################################################################################################### 
########### loading corrected p-values, output of the script "coexpression_testing_randomization" #########
###########################################################################################################

load(verbose=T,file=paste0("intermediate_files/peprec.pval.randomisedcorrection.mat.",pdum.or.mouse,".RData"))

limit.pval.used <- limit.pval.corrected # 0.05

###########################################################################################
########### we define here the significant pairs of markers which are coexpressed #########
###########################################################################################

signif.coexpr <- lapply(1:nrow(pval.used), function(x) {l<-which(pval.used[x,]<limit.pval.used);return(l[which(l>x)])} )
sign.peps <- lapply(1:nrow(pval.used), function(x) peprec.pep.mat[x, signif.coexpr[[x]]])
sign.recs <- lapply(1:nrow(pval.used), function(x) peprec.rec.mat[x, signif.coexpr[[x]]])
pep.to.be.considered.incoexpr <- which(unlist(lapply(1:length(sign.peps), function(x) length(sign.peps[[x]])!=0)))
list.signif <- lapply(pep.to.be.considered.incoexpr, function(x) cbind(peprec.pep.mat.commonname[x, signif.coexpr[[x]]], peprec.rec.mat.commonname[x, signif.coexpr[[x]]], peprec.corr.mat[x, signif.coexpr[[x]]], pval.used[x, signif.coexpr[[x]]], peprec.or.mat[x, signif.coexpr[[x]]], peprec.pval.fisher.mat[x, signif.coexpr[[x]]], peprec.numintersect.mat[x, signif.coexpr[[x]]], peprec.numcellnon0.1.mat[x, signif.coexpr[[x]]],peprec.numcellnon0.2.mat[x, signif.coexpr[[x]]], ncol(tab.current), peprec.mic.mat[x, signif.coexpr[[x]]]))
coexpr.vec <- do.call(rbind, list.signif)
colnames(coexpr.vec) <- c("name.1", "name.2", "correlation", "pval correlation", "OR Fisher", "Pval Fisher", "num.cells.coexpress", "num.cells.group1", "num.cells.group2","total.num.cells", "MIC")

################################################################################################
########### outputs the results - we save results.analyses to use later with corrected pvals ####
################################################################################################

write.table(coexpr.vec, paste0("intermediate_files/allmarkersvsallmarkers_",pdum.or.mouse,"_",type.of.mouse.data,"_coexpression_",randomisation_folder.mouse,".txt"), sep="\t", col.names = NA)
results.analyses <- list(signif.coexpr, peprec.pep.mat.commonname, peprec.rec.mat.commonname, peprec.corr.mat, peprec.pval.mat, peprec.or.mat, peprec.pval.fisher.mat, peprec.numintersect.mat, peprec.numcellnon0.1.mat,peprec.numcellnon0.2.mat, ncol(tab.current), peprec.mic.mat)
save(results.analyses, file=paste0("intermediate_files/allmarkersvsallmarkers_",pdum.or.mouse,"_",type.of.mouse.data,"_coexpression_",randomisation_folder.mouse,".RData") )


################################################################################################
########### fixing the correlation matrix - diagonal and same vs. same markers are set to 1 and 
########### lines and columns (excluding the diagonal) from markers which are not expressed are 
########### set to 0  ##########################################################################
################################################################################################

fix_corr_matrix <- function(corr.mat){
  
  ind.genes.notfound <- which(unlist(lapply(1:nrow(corr.mat), function(x) corr.mat[x,x]<0.99)))
  res <- corr.mat
  for (i in 1:length(ind.genes.notfound)){
    for (j in 1:ncol(corr.mat)){
      res[ind.genes.notfound[i], j] <- 0
    }
    res[ind.genes.notfound[i], ind.genes.notfound[i]] <- 1
  }
  
  for (j in 1:length(ind.genes.notfound)){
    for (i in 1:nrow(corr.mat)){
      res[i, ind.genes.notfound[j]] <- 0
    }
    res[ind.genes.notfound[j], ind.genes.notfound[j]] <- 1
  }
  return(res)
}
peprec.corr.mat <- fix_corr_matrix(peprec.corr.mat)

################################################################################################
########### saving the correlation matrix ###############################################
################################################################################################

save(peprec.corr.mat, file=paste0("intermediate_files/peprec.corr.mat.",pdum.or.mouse,".RData"))

##############################################################################################
################# what follows are analyses based on saved files for correction of pvalues ###
##############################################################################################

########################################
############### parameters #############
########################################

pdum.or.mouse<-"mouse"  # c("pdum", "mouse")

########################################################################################
############### loading both the mouse and platynereis correlation matrices ############
########################################################################################

# platynereis
load("intermediate_files/peprec.corr.mat.pdum.RData", verbose=T)
peprec.corr.mat <- fix_corr_matrix(peprec.corr.mat)
corr.vectorised <- do.call(c, lapply(1:(ncol(peprec.corr.mat)-1), function(x) peprec.corr.mat[(x+1):nrow(peprec.corr.mat),x]))
corr.vectorised.pdum <- corr.vectorised
# mouse
load("intermediate_files/peprec.corr.mat.mouse.RData", verbose=T)
peprec.corr.mat <- fix_corr_matrix(peprec.corr.mat)
corr.vectorised <- do.call(c, lapply(1:(ncol(peprec.corr.mat)-1), function(x) peprec.corr.mat[(x+1):nrow(peprec.corr.mat),x]))
corr.vectorised.mouse <- corr.vectorised

# this tests the global correlation between the two species (on the upper left triangular matrix to exclude same vs same data points)
cor.test(corr.vectorised.pdum, corr.vectorised.mouse) ### the "correlation_pdum_mouse_romanov.txt" file

##########################################################################################################################
############### loading both the mouse and platynereis results and corrected pvalues (output of coexpression_testing_randomization) ############
##########################################################################################################################

load(verbose=T, paste0("intermediate_files/peprec.pval.randomisedcorrection.mat.pdum_.RData"))
load(verbose=T, paste0("intermediate_files/allmarkersvsallmarkers_pdum_coexpression.RData") )
results.analyses.1 <- results.analyses
peprec.pep.mat.commonname.1 <- results.analyses[[2]]
peprec.rec.mat.commonname.1 <- results.analyses[[3]] # this is just the transpose of peprec.pep.mat.commonname.1
peprec.corr.mat.1 <- results.analyses[[4]]
peprec.pval.mat.1 <- results.analyses[[5]]
peprec.correctedpval.mat.1 <- pval.randomised.correction
peprec.pval.fisher.mat.1 <- results.analyses[[7]]
load(verbose=T, paste0("intermediate_files/peprec.pval.randomisedcorrection.mat.mouse.RData"))
load(verbose=T, paste0("intermediate_files/allmarkersvsallmarkers_mouse_coexpression.RData") )
results.analyses.2 <- results.analyses
peprec.pep.mat.commonname.2 <- results.analyses[[2]]
peprec.rec.mat.commonname.2 <- results.analyses[[3]] # this is just the transpose of peprec.pep.mat.commonname.2
peprec.corr.mat.2 <- results.analyses[[4]]
peprec.pval.mat.2 <- results.analyses[[5]]
peprec.correctedpval.mat.2 <- pval.randomised.correction
peprec.pval.fisher.mat.2 <- results.analyses[[7]]


peprec.corr.mat.1 <- fix_corr_matrix(peprec.corr.mat.1)
peprec.corr.mat.2 <- fix_corr_matrix(peprec.corr.mat.2)
for (i in 1:nrow(peprec.pval.mat.1)){
  peprec.pval.mat.1[i,i] <- 0
  peprec.pval.mat.2[i,i] <- 0
}

#############################################################################
############### here we define the coexpressed pairs that are conserved #####
#############################################################################

inds.conserved.corrs <- lapply(1:(ncol(peprec.corr.mat.1)-1), function(x) which(peprec.correctedpval.mat.1[x, (x+1):ncol(peprec.correctedpval.mat.1)] < limit.pval.used & peprec.correctedpval.mat.2[x, (x+1):ncol(peprec.correctedpval.mat.2)] < limit.pval.used)+x)
pairs.coexpr <- list()
t<-1
for (i in 1:length(inds.conserved.corrs)){
  if (length(inds.conserved.corrs[[i]]>0)){
    for (j in 1:length(inds.conserved.corrs[[i]])){
      pairs.coexpr[[t]] <- c(i, inds.conserved.corrs[[i]][j])
      t<-t+1
    }
  }
}

output.signif <- do.call(rbind,lapply(1:length(pairs.coexpr), function(x) cbind(peprec.pep.mat.commonname.1[pairs.coexpr[[x]][1],1], peprec.pep.mat.commonname.1[pairs.coexpr[[x]][2],1], peprec.corr.mat.1[pairs.coexpr[[x]][1],pairs.coexpr[[x]][2]],  peprec.pval.mat.1[pairs.coexpr[[x]][1],pairs.coexpr[[x]][2]], peprec.correctedpval.mat.1[pairs.coexpr[[x]][1],pairs.coexpr[[x]][2]], peprec.pval.fisher.mat.1[pairs.coexpr[[x]][1],pairs.coexpr[[x]][2]], peprec.pep.mat.commonname.2[pairs.coexpr[[x]][1],pairs.coexpr[[x]][2]], peprec.pep.mat.commonname.2[pairs.coexpr[[x]][2],1], peprec.corr.mat.2[pairs.coexpr[[x]][1],pairs.coexpr[[x]][2]], peprec.pval.mat.2[pairs.coexpr[[x]][1],pairs.coexpr[[x]][2]], peprec.correctedpval.mat.2[pairs.coexpr[[x]][1],pairs.coexpr[[x]][2]], peprec.pval.fisher.mat.2[pairs.coexpr[[x]][1],pairs.coexpr[[x]][2]])) )
output.signif.df <- data.frame(output.signif)
colnames(output.signif.df) <- c("name.Pdum.1","name.Pdum.2","corr.Pdum", "Pval.Pdum", "corrected_Pval.Pdum", "FTest_Pval.Pdum", "name.mouse.1","name.mouse.2","corr.mouse", "Pval.mouse", "corrected_Pval.mouse", "FTest_Pval.mouse")

##################################################################################
############### outputs the final table of conserved markers #####################
##################################################################################

write.table(col.names=NA, output.signif.df, paste0("figures/conserved_coexpression_pdum_romanov.txt"), sep="\t")

##################################################################################
############### the supplemental figures from the Williams et al. 2017 paper #####
##################################################################################

rownames(peprec.corr.mat.1) <- peprec.rec.mat.commonname.1[1,]
colnames(peprec.corr.mat.1) <- peprec.rec.mat.commonname.1[1,]
par(cex=1)
pdf(paste0("figures/correlations_pdum_allmarkers.pdf"))
corrplot(fix.corr.matrix(peprec.corr.mat.1)[1:28, ], p.mat = peprec.pval.mat.1[1:28, ], sig.level = limit.pval.used, insig = "blank", method = "circle", type = "upper", tl.cex=0.3)
dev.off()

rownames(peprec.corr.mat.2) <- peprec.rec.mat.commonname.2[1,]
colnames(peprec.corr.mat.2) <- peprec.rec.mat.commonname.2[1,]
par(cex=1)
pdf(paste0("figures/correlations_mouse_allmarkers.pdf"))
corrplot(fix.corr.matrix(peprec.corr.mat.2)[1:28, ], p.mat = peprec.pval.mat.2[1:28, ], sig.level = limit.pval.used, insig = "blank", method = "circle", type = "upper", tl.cex=0.3)
dev.off()

