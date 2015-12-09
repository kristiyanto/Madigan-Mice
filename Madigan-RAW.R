# This is for class TCSS 592 
# By Daniel Kristiyanto (danielkr@uw.edu)
# Fall 2015
# LATEST VERSION.
# THIS FILE READ THE CEL FILE PROVIDED BY MADIGAN 
# (RAW DATA)
#################################### REQUIRED PACKAGES ######################################

# Setting environtment
source("https://bioconductor.org/biocLite.R")
library("oligo")
library("limma")
library("gplots")
library("reshape2")
library("RColorBrewer")
library("fields")
library("mogene21sttranscriptcluster.db")
working.dir <- "~/Google Drive/BIOINFORMATICS/Madigan"
setwd(working.dir)

#################################### PARAMETERS ######################################
# Export to tableu
tableu      <- 0
graphs      <- 0
readfromraw <- 0
save.this.session <- 0
outfile.id  <- "CEL file"

group1 <- rep(13:24)
group2 <- rep(1:12)
group3 <- rep(25:36)
group4 <- rep(37:48)
group5 <- rep(49:60)
group6 <- rep(61:72)
group7 <- rep(73:84)
group8 <- rep(85:96)

group.names       <- c("S","S","M","M","B","B","MB","MB")
subgroup.names    <- rep(c(".P",".I"),4)
all.groups        <- NULL
all.grp           <- NULL
for(n in 1:8)
{
  x <- paste0(group.names[n],subgroup.names[n])
  x <- rep(x,13)
  x <- make.unique(x)
  all.groups <- c(all.groups,x[2:13])
  x <- paste0("GROUP",n)
  x <- rep(x,13)
  x <- make.unique(x)
  all.grp <- c(all.grp,x[2:13])
}
###################################### RAW DATA ######################################
if(readfromraw==1)
{
  # Read the files
  cel.files     <- list.celfiles(path="~/Google Drive/BIOINFORMATICS/Madigan/LPS Study",full.names=TRUE)
  core.cel      <- read.celfiles(cel.files)
  core.qsp <- normalize(core.cel,method="qspline", copy=TRUE, verbose=T, target='core')
  core.rma <- rma(core.qsp, background=T, normalize=F, subset=NULL, target="probeset")
  
  # Read Annotation Data
  # Compare Genes 3 & 4
  mogene21sttranscriptcluster()
  ann.data   <- data.frame(ENTREZ=sapply(contents(mogene21sttranscriptclusterENTREZID),paste,collapse=", "),
                           SYMBOL=sapply(contents(mogene21sttranscriptclusterSYMBOL), paste, collapse=", "), 
                           DESC=sapply(contents(mogene21sttranscriptclusterGENENAME), paste, collapse=", "))
  # core.data.ann       <- merge(core.data.norm, annotation.mogene, by.x=0, by.y=0, all=T)
  save.image("celnormalized.RData")
}else{
  load("celnormalized.RData")
}

#################################### DIFFERENTIAL EXPRESSION ######################################
    
    p                     <- 0.05
    curr.diff             <- exprs(core.rma)
    colnames(curr.diff)   <- all.grp
    
    D <- as.factor(c(rep("GROUP1",12),rep("GROUP2",12),rep("GROUP3",12),rep("GROUP4",12),rep("GROUP5",12),
           rep("GROUP6",12),rep("GROUP7",12),rep("GROUP8",12)))
    
    design                <- model.matrix(~0 + D)
    colnames(design)      <- levels(D)
    
    fit                   <- lmFit(curr.diff,design) 
    cont.matrix           <- makeContrasts(GROUP1-GROUP2,GROUP3-GROUP4,GROUP5-GROUP6,GROUP7-GROUP8, levels=design)
    fit2                  <- contrasts.fit(fit,cont.matrix)
    fit2                  <- eBayes(fit2)
    topTable(fit2,p.value=p, adjust="fdr")
    
    results              <- decideTests(fit2,adjust.method="fdr",p.value=p)
    summary(results)
    selected             <- p.adjust(fit2$F.p.value, method = "fdr") < p  # P Adj ####
    sum(selected)
    
    side.col.full   <- NULL
    for(c in brewer.pal(8,"Spectral"))
      side.col.full <- c(side.col.full,rep(c,12))

    heatmap.2(as.matrix(curr.diff[selected==1,]), main = "GLOBAL. F.p.value<5%, Bon", key=TRUE, symkey=FALSE, density.info="none", trace="none",
                ColSideColors = side.col.full,margins =c(12,9),
                col=rev(brewer.pal(11, "Spectral")), scale = "none")

for(h in c(1,2,3,4))
{
  if(h==1){s = c(seq(1,12),seq(13,24))} else if (h==2) {s = seq(25,48)} else if (h==3) {s = seq(49,72)} else {s=seq(73:96)}
  side.col.full <- c(rep("red",12),rep("green",12))
  
  cont = as.factor(c("GROUP1-GROUP2","GROUP3-GROUP4","GROUP5-GROUP6","GROUP7-GROUP8"))
  heatmap.2(as.matrix(curr.diff[selected,s]), main = cont[h], key=TRUE, symkey=FALSE, density.info="none", trace="none",
         ColSideColors = side.col.full,margins =c(12,9),
         col=rev(brewer.pal(11, "Spectral")), scale = "none")

  sg             <- p.adjust(fit2$p.value[,h], method = "fdr") < p
  if(sum(sg)>2)
  heatmap.2(as.matrix(curr.diff[sg,s]), main = paste0(cont[h],"Locally Diff Exp"), key=TRUE, symkey=FALSE, density.info="none", trace="none",
          ColSideColors = side.col.full,margins =c(12,9),
          col=rev(brewer.pal(11, "Spectral")), scale = "none")
}

if(save.this.session==1)
save.image("tomarkdown.RData")
    