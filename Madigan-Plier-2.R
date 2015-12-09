# This is for class TCSS 592 
# By Daniel Kristiyanto (danielkr@uw.edu)
# Fall 2015
# LATEST VERSION.
# THIS FILE READ THE CSV FILE PROVIDED BY MADIGAN 
# (PLIER NORMALIZED)
#################################### REQUIRED PACKAGES ######################################

# Setting environtment
source("https://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("affxparser")
# library("affy")
library("affxparser")
library("oligo")
library("limma")
library("gplots")
library("scales")
library("reshape2")
working.dir <- "~/Google Drive/BIOINFORMATICS/Madigan"
setwd(working.dir)

#################################### PARAMETERS ######################################
# Export to tableu
tableu      <- 0
graphs      <- 0
readfromraw <- 0
outfile.id  <- "gene-level"

group1 <- rep(13:24)
group2 <- rep(1:12)
group3 <- rep(25:36)
group4 <- rep(37:48)
group5 <- rep(49:60)
group6 <- rep(61:72)
group7 <- rep(73:84)
group8 <- rep(85:96)

group.names       <- c("S","S","M","M","B","B","MB","MB")
subgroup.names    <- rep(c("-P","-I"),4)
all.groups        <- NULL

for(n in 1:8)
{
  x <- paste0(group.names[n],subgroup.names[n])
  x <- rep(x,13)
  x <- make.unique(x)
  all.groups <- c(all.groups,x[2:13])
}

###################################### READING THE FILES (FROM KEVIN) ######################################
if(readfromraw==1)
{
  # Read the files
  col.names       <- NULL
  row.names       <- NULL
  core.data        <- NULL
  
  setwd(paste0(working.dir,"/PLIER/Madigan_LPS_CHP/"))
  ## LOOP THROUGH ALL FILES AND APPEND THE COLUMNS
  for (file in list.files())
  {
    cat("Working on", file,"\n")
    shortname   <- sub("\\..*","", file)
    col.names   <- c(col.names, shortname)
    curr.chp    <- readChp(file, withQuant = TRUE)
    
    if (length(row.names) == 0) # THIS IS THE FIRST TIME SO CAPTURE ROW NAMES
    row.names   <- curr.chp$QuantificationEntries[[1]] 
    core.data   <- cbind(core.data, curr.chp$QuantificationEntries[[2]])
  }
  
  core.data           <- as.data.frame(core.data)
  rownames(core.data) <- row.names
  colnames(core.data) <- all.groups
  
  core.data.norm      <- normalizeCyclicLoess(as.matrix(core.data[,-13]),method="affy", iterations=2)
  core.data.norm      <- cbind(core.data.norm, core.data[,13])
  colnames(core.data.norm)[96] <- "S-I.1"
  core.data.norm      <- core.data.norm[,all.groups]
  write.csv(core.data.norm,row.names = T, file="LOESS-NORMALIZED.csv")
  save.image("normalized.RData")
}else{
  load("normalized.RData")
}
#################################### ANNOTATION ######################################
# Annotation Data
# Compare Genes 3 & 4
# biocLite("mogene21sttranscriptcluster.db")
library("mogene21sttranscriptcluster.db")
mogene21sttranscriptcluster()
annotation.mogene   <- data.frame(ENTREZ=sapply(contents(mogene21sttranscriptclusterENTREZID),paste,collapse=", "),
                    ACCNUM=sapply(contents(mogene21sttranscriptclusterACCNUM), paste, collapse=", "), 
                    SYMBOL=sapply(contents(mogene21sttranscriptclusterSYMBOL), paste, collapse=", "), 
                    DESC=sapply(contents(mogene21sttranscriptclusterGENENAME), paste, collapse=", "))
# write.csv(core.data,file="OUTPUT/core.data.txt",sep="\t", row.names = T)
core.data.ann       <- merge(core.data.norm, annotation.mogene, by.x=0, by.y=0, all=T)
#################################### DIFFERENTIAL EXPRESSION ######################################

hmap          <- 1

if(hmap==1)
{
  
  p.distribution            <- NULL
  t.distribution            <- NULL
  selected.genes            <- NULL
  for(n in c(1,3,5,7))
  {
      
      a                     <- eval(parse(text=paste0("group",n)))
      b                     <- eval(parse(text=paste0("group",n+1)))
      grp.a                 <- paste0("group",n)
      grp.b                 <- paste0("group",n+1)
      graph.name            <- paste(grp.a,grp.b,"CHP") 
      p                     <- 0.05
      
      curr.diff <- as.data.frame(core.data.ann[core.data.ann$SYMBOL!="NA",c(a,b)])
      curr.ori   <- as.data.frame(core.data[core.data.ann$SYMBOL!="NA",c(a,b)])
      if(n==1)
      {
        pert.names          <- factor(c(rep(paste0("group",n),12),rep(paste0("group",n+1),11)))
        curr.diff[,13]      <- NULL
        curr.ori[,13]       <- NULL
      } else {
        pert.names          <- factor(c(rep(paste0("group",n),12),rep(paste0("group",n+1),12)))
      }
      
      design                <- model.matrix(~0 + pert.names)
      colnames(design)      <- levels(pert.names)
    
      fit                   <- lmFit(curr.diff,design) 
      cont.matrix           <- makeContrasts(eval(parse(text=grp.a))-eval(parse(text=grp.b)), levels=design)
      fit2                  <- contrasts.fit(fit,cont.matrix)
      fit2                  <- eBayes(fit2)
      topTable(fit2,p.value=0.05, adjust="fdr")
      core.annotated[row.names(topTable(fit2,p.value = 0.01, adjust="bonferroni")),]
      
      selected              <- p.adjust(fit2$p.value, method = "bonferroni") < p  # P ADJ ####
      
      selected.genes        <- rbind(selected.genes, cbind(graph.name,
                                                           core.data.ann[row.names(core.data.ann) %in% row.names(curr.diff[selected==1,]),]))
      row_distance          <- dist(as.matrix(curr.diff[selected==1,]), method = "manhattan")
      row_cluster           <- hclust(row_distance, method = "ward.D")
      col_distance          <- dist(t(as.matrix(curr.diff[selected==1,])), method = "manhattan")
      col_cluster           <- hclust(col_distance, method = "ward.D")
      
      
      png(file=paste0(working.dir,"/graph/",outfile.id,".",graph.name,".volcano.png"),width = 900, height = 900, units = "px")
      volcanoplot(fit2, main=graph.name, highlight = 100)
      dev.off()
      
      library("RColorBrewer")
      png(file=paste0(working.dir,"/graph/",outfile.id,".",graph.name,".heatmap.png"),width = 900, height = 1200, units = "px")
      if(n==1)
      {
        heatmap.2(as.matrix(curr.diff[selected < 0.95,]), main = graph.name, key=TRUE, symkey=FALSE, density.info="none", trace="none", 
                  Rowv = as.dendrogram(row_cluster), ColSideColors = c(rep("green", 12),rep("blue", 11)),margins =c(12,9),
                  Colv = as.dendrogram(col_cluster),col=brewer.pal(6, "Spectral"), scale = "row")
      }else {
        heatmap.2(as.matrix(curr.diff[selected==1,]), main = graph.name, key=TRUE, symkey=FALSE, density.info="none", trace="none",
                  Rowv = as.dendrogram(row_cluster), ColSideColors = c(rep("green", 12),rep("blue", 12)),margins =c(12,9),
                  Colv = as.dendrogram(col_cluster),col=rev(brewer.pal(11, "Spectral")), scale = "row")
      }
      dev.off()
      
      
      png(file=paste0(working.dir,"/graph/",outfile.id,".",graph.name,".p.png"),width = 900, height = 900, units = "px")
      hist(fit2$p.value, main=paste("P Value",graph.name), freq=T, breaks=seq(0,1,p), xlab="PValue", col=brewer.pal(11,"Spectral"))
      dev.off()
      
      png(file=paste0(working.dir,"/graph/",outfile.id,".",graph.name,".t.png"),width = 900, height = 900, units = "px")
      hist(fit2$t, main=paste("T Statistic",graph.name), freq=T, breaks=seq(-10,10,0.5), xlab="TStatistic", col=brewer.pal(12, "Paired"))
      dev.off()
    p.distribution        <- rbind(p.distribution, as.data.frame(cbind("p.value", graph.name,row.names(fit2$p.value),fit2$p.value)))
    t.distribution        <- rbind(t.distribution, as.data.frame(cbind("t-test",graph.name,row.names(fit2$t),fit2$t)))  
  }
  p.tableu          <- as.data.frame(rbind(p.distribution,t.distribution))
  names(p.tableu)   <- c("LABEL","Group","Probe No","Value")
  write.csv(selected.genes,file=paste0(working.dir,"/OUTPUT/",outfile.id,".selected-genes.csv"))
  
  
  #################################### SELECTED GENES ######################################
  SG <- selected.genes
  gr12 <- SG[SG$graph.name==levels(SG$graph.name)[1], "SYMBOL"]
  gr34 <- SG[SG$graph.name==levels(SG$graph.name)[2], "SYMBOL"]
  gr56 <- SG[SG$graph.name==levels(SG$graph.name)[3], "SYMBOL"]
  gr78 <- SG[SG$graph.name==levels(SG$graph.name)[4], "SYMBOL"]
  
  aha.genes <- Reduce(intersect, list(gr12,gr34,gr56,gr78))
  
  small.genes <- (SG[SG$SYMBOL %in% aha.genes,])
  
  
  for(h in 1:4)
  {
  
    if(h==1){s = c(seq(1,12),seq(13,24))} else if (h==2) {s = seq(25,48)} else if (h==3) {s = seq(49,72)} else {s=seq(73:96)}
  
  g           <- levels(small.genes$graph.name)[h]
  this.set    <- small.genes[small.genes$graph.name==g,c(all.groups,"SYMBOL")]
  row.names(this.set) <- make.unique(as.character(this.set$SYMBOL))
  this.set    <- this.set[,all.groups]
  this.set    <- this.set[,s]
  png(file=paste0(working.dir,"/graph/4sel-",g,".png"),width = 900, height = 900, units = "px")
  heatmap.2(as.matrix(this.set), main = g, key=TRUE, symkey=FALSE, density.info="none", trace="none",
            margins =c(12,9), dendrogram = "column", ColSideColors = c(rep("green", 12),rep("blue", ncol(this.set)-12))
            ,col=rev(brewer.pal(11, "Spectral")),scale="row")
  dev.off()
  }
}

#################################### ALL GROUPS ######################################
all.data =1
if(all.data==1)
{
  
    p                     <- 0.05

    curr.diff  <- as.data.frame(core.data.ann[,all.groups[-13]])

    D <- as.factor(c(rep("GROUP1",12),rep("GROUP2",11),rep("GROUP3",12),rep("GROUP4",12),rep("GROUP5",12),
           rep("GROUP6",12),rep("GROUP7",12),rep("GROUP8",12)))
    
    design                <- model.matrix(~0 + D)
    colnames(design)      <- levels(D)
    
    fit                   <- lmFit(curr.diff,design) 
    cont.matrix           <- makeContrasts(GROUP1-GROUP2,GROUP3-GROUP4,GROUP5-GROUP6,GROUP7-GROUP8, levels=design)
    fit2                  <- contrasts.fit(fit,cont.matrix)
    fit2                  <- eBayes(fit2)
    topTable(fit2,coef=1, adjust="bonferroni")
    
    # results              <- decideTests(fit2,method="separate",adjust.method="bonferroni",p.value=p)

    selected              <- p.adjust(fit2$F.p.value, method = "bonferroni") < p
    row_distance          <- dist(as.matrix(curr.diff[selected==1,]), method = "manhattan")
    row_cluster           <- hclust(row_distance, method = "ward.D")
    col_distance          <- dist(t(as.matrix(curr.diff[selected==1,])), method = "manhattan")
    col_cluster           <- hclust(col_distance, method = "ward.D")
    
    
    png(file=paste0(working.dir,"/graph/",outfile.id,".",graph.name,".volcano.png"),width = 900, height = 900, units = "px")
    volcanoplot(fit2, highlight = 100)
    dev.off()
    
    library("RColorBrewer")
    side.col.full <- NULL
    for(c in brewer.pal(8,"Spectral"))
      side.col.full <- c(side.col.full,rep(c,12))
      side.col.full <- side.col.full[-13]
    
    png(file=paste0(working.dir,"/graph/",outfile.id,".globalheat.png"),width = 900, height = 1200, units = "px")
      heatmap.2(as.matrix(curr.diff[selected==1,]), main = "GLOBAL. F.p.value<5%, Bon", key=TRUE, symkey=FALSE, density.info="none", trace="none",
                Rowv = as.dendrogram(row_cluster), ColSideColors = side.col.full,margins =c(12,9),
                Colv = as.dendrogram(col_cluster),col=rev(brewer.pal(11, "Spectral")), scale = "row")
    dev.off()
    
    results <- decideTests(fit2, adjust.method = "bonferroni")
    
    png(file=paste0(working.dir,"/graph/",outfile.id,"global.Fp.png"),width = 900, height = 900, units = "px")
    hist(fit2$F.p.value, main=paste("All. F. P Value"), freq=T, breaks=seq(0,1,p), xlab="PValue", col=brewer.pal(11,"Spectral"))
    dev.off()
    
    png(file=paste0(working.dir,"/graph/",outfile.id,"Global.t.png"),width = 900, height = 900, units = "px")
    hist(fit2$t, main=paste("All T Statistic"), freq=T, breaks=seq(-10,10,0.5), xlab="TStatistic", col=brewer.pal(12, "Paired"))
    dev.off()

    for(h in c(1,3,4))
    {
      
      if(h==1){s = c(seq(1,12),seq(13,24))} else if (h==2) {s = seq(25,48)} else if (h==3) {s = seq(49,72)} else {s=seq(73:96)}
      side.col.full <- c(rep("red",12),rep("green",length(s)-12))
      
    png(file=paste0(working.dir,"/graph/",outfile.id,h,".globalheat.png"),width = 900, height = 1200, units = "px")
    heatmap.2(as.matrix(), main = h, key=TRUE, symkey=FALSE, density.info="none", trace="none",
             ColSideColors = side.col.full,margins =c(12,9),
             col=rev(brewer.pal(11, "Spectral")), scale = "none")
    dev.off()
    }

}


#################################### EXPORT TO TABLEU ######################################

# TABLEU #########
# Augmenting the data for Tableu visualization 
if(tableu==1)
{
  cn       <- merge(core.data.norm, annotation.mogene, by.x=0, by.y=0, all=T)
  core.annotated      <- NULL
  for(i in 1:8)
  {
    SQ          <- group.names[i]
    IU          <- subgroup.names[i]
    assign(paste(text=paste0("subset",i)), cbind(paste0("GROUP",i),SQ,IU,core.data.ann[core.data.ann$SYMBOL!="NA",c("ENTREZ","ACCNUM","SYMBOL","DESC")],
                                                 core.data.ann[core.data.ann$SYMBOL!="NA",eval(parse(text=paste0("group",i)))]))
    
    tmp         <- get(paste(text=paste0("subset",i)))
    names(tmp)  <- c("GROUP","SQ","IU","ENTREZ","ACCNUM","SYMBOL","DESC",
                     "MOUSE1","MOUSE2","MOUSE3","MOUSE4","MOUSE5","MOUSE6","MOUSE7",
                     "MOUSE8","MOUSE9","MOUSE10","MOUSE11","MOUSE12")
    core.annotated <- rbind(core.annotated, tmp)
    if(i==2)
    {
      core.annotated$MOUSE1 <- NA
    }
    
  }
  
  dim(core.annotated)
  core.reshape  <- melt(core.annotated, id.vars = c("GROUP","ENTREZ","SQ","IU","ACCNUM","SYMBOL","DESC"),na.rm = T)
  colnames(core.reshape)[8] <- "EXPERIMENT"
  colnames(core.reshape)[9] <- "EXPRESSION"
  write.csv(core.reshape,file=paste0(working.dir,"/OUTPUT/",outfile.id,"core.txt"), row.names = T)
  
  
  ### SELECTED GENES ##
  selected.reshape <- melt(selected.genes, id.vars = c("graph.name", "Row.names","ENTREZ","ACCNUM","SYMBOL","DESC"))
  write.csv(selected.reshape,file=paste0(working.dir,"/OUTPUT/",outfile.id,"selected.tableu.txt"), row.names = T)
}
