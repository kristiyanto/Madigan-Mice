# This is for class TCSS 592 
# By Daniel Kristiyanto (danielkr@uw.edu)
# Fall 2015
#################################### REQUIRED PACKAGES ######################################

working.dir <- "~/Google Drive/BIOINFORMATICS/Madigan"
setwd(working.dir)

#################################### PARAMETERS ######################################
# Export to tableu
tableu <- 1
graphs <- 0
diffex.type <- 1


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

###################################### RAW DATA ######################################
# THIS SECTION IS TO NORMALIZE THE RAW DATA AND SAVE IT      #
#

# Setting environtment
source("https://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("affxparser")
library("affy")
library("affxparser")

library("oligo")
library("limma")
library("gplots")
library("scales")
library("reshape2")

###################################### READING THE FILES (FROM KEVIN) ######################################

# Read the files
col.names       <- NULL
row.names       <- NULL
#READ A FILE 

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

#################################### ANNOTATION ######################################
# Annotation Data
# Compare Genes 3 & 4
# biocLite("mogene21sttranscriptcluster.db")
# library("mogene21sttranscriptcluster.db")
mogene21sttranscriptcluster()
annotation.mogene   <- data.frame(ENTREZ=sapply(contents(mogene21sttranscriptclusterENTREZID),paste,collapse=", "),
                    ACCNUM=sapply(contents(mogene21sttranscriptclusterACCNUM), paste, collapse=", "), 
                    SYMBOL=sapply(contents(mogene21sttranscriptclusterSYMBOL), paste, collapse=", "), 
                    DESC=sapply(contents(mogene21sttranscriptclusterGENENAME), paste, collapse=", "))
# write.csv(core.data,file="OUTPUT/core.data.txt",sep="\t", row.names = T)
core.data.ann       <- merge(core.data, annotation.mogene, by.x=0, by.y=0, all=T)

#################################### DIFFERENTIAL EXPRESSION ######################################

hmap          <- 1

if(hmap==1)
{
  table.for.diff            <- core.data
  p.distribution            <- NULL
  t.distribution            <- NULL
  selected.genes            <- NULL
  
  
  for(n in c(1,3,5,7))
  {
     
      a                     <- eval(parse(text=paste0("group",n)))
      b                     <- eval(parse(text=paste0("group",n+1)))
      curr.diff             <- table.for.diff[,c(a,b)]
      grp.a                 <- paste0("group",n)
      grp.b                 <- paste0("group",n+1)
      
      pert.names            <- factor(c(rep(paste0("group",n),12),rep(paste0("group",n+1),12)))

      design                <- model.matrix(~0 + pert.names)
      colnames(design)      <- levels(pert.names)
      graph.name            <- paste(grp.a,grp.b,"PLIER") 
    
      fit                   <- lmFit(curr.diff,design) 
      cont.matrix           <- makeContrasts(eval(parse(text=grp.a))-eval(parse(text=grp.b)), levels=design)
      fit2                  <- contrasts.fit(fit,cont.matrix)
      fit2                  <- eBayes(fit2)
      topTable(fit2,coef=1, adjust="bonferroni")
      
      results               <- decideTests(fit2)
      
      png(file=paste0(working.dir,"/graph/venn-",graph.name,".png"),width = 900, height = 900, units = "px")
      vennDiagram(results, main=graph.name, names="")
      dev.off()
      
      png(file=paste0(working.dir,"/graph/volcano-",graph.name,".png"),width = 900, height = 900, units = "px")
      volcanoplot(fit2, main=graph.name, names=fit2$genes$ID, highlight = 50)
      dev.off()
      
      selected              <- p.adjust(fit2$p.value < 0.01, method = "bonferroni")
      selected.genes        <- rbind(selected.genes, cbind(graph.name,core.data.ann[selected==1,]))
      row_distance          <- dist(as.matrix(curr.diff[selected==1,]), method = "manhattan")
      row_cluster           <- hclust(row_distance, method = "ward.D")
      col_distance          <- dist(t(as.matrix(curr.diff[selected==1,])), method = "manhattan")
      col_cluster           <- hclust(col_distance, method = "ward.D")
      
      png(file=paste0(working.dir,"/graph/",graph.name,".png"),width = 900, height = 900, units = "px")
      heatmap.2(as.matrix(curr.diff[selected==1,]), main = graph.name, key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="none",
                Rowv = as.dendrogram(row_cluster), ColSideColors = c(rep("green", 12),rep("blue", 12)),margins =c(12,9),
                Colv = as.dendrogram(col_cluster))
      dev.off()
      
      p.distribution        <- rbind(p.distribution, as.data.frame(cbind("p.value", graph.name,row.names(fit2$p.value),fit2$p.value)))
      t.distribution        <- rbind(t.distribution, as.data.frame(cbind("t-test",graph.name,row.names(fit2$t),fit2$t)))  
  }
  p.tableu          <- as.data.frame(rbind(p.distribution,t.distribution))
  names(p.tableu)   <- c("LABEL","Group","Probe No","Value")
  write.csv(p.tableu,file=paste0(working.dir,"/OUTPUT/p-dist-plier.csv"), row.names = F)
  write.csv(selected.genes,file=paste0(working.dir,"/OUTPUT/selectedgenes-plier.csv"))
}



# TABLEU #########
# Augmenting the data for Tableu visualization 
tableu = 1
if(tableu==1)
{
  
  core.annotated      <- NULL
  for(i in 1:8)
  {
    SQ          <- group.names[i]
    IU          <- subgroup.names[i]
    assign(paste(text=paste0("subset",i)), cbind(paste0("GROUP",i),SQ,IU,core.data.ann[,c("ENTREZ","ACCNUM","SYMBOL","DESC")],
                                                 core.data.ann[,eval(parse(text=paste0("group",i)))]))
    
    tmp         <- get(paste(text=paste0("subset",i)))
    names(tmp)  <- c("GROUP","SQ","IU","ENTREZ","ACCNUM","SYMBOL","DESC",
                     "MOUSE1","MOUSE2","MOUSE3","MOUSE4","MOUSE5","MOUSE6","MOUSE7",
                     "MOUSE8","MOUSE9","MOUSE10","MOUSE11","MOUSE12")
    core.annotated <- rbind(core.annotated, tmp)
  }
  
  dim(core.annotated)
  core.reshape  <- melt(core.annotated, id.vars = c("GROUP","ENTREZ","SQ","IU","ACCNUM","SYMBOL","DESC"))
  colnames(core.reshape)[8] <- "EXPERIMENT"
  colnames(core.reshape)[9] <- "EXPRESSION"
  write.csv(core.reshape,file=paste0(working.dir,"/OUTPUT/core.aug.plier.txt"), row.names = T)
  
  
  ### SELECTED GENES ##
  selected.reshape <- melt(selected.genes, id.vars = c("graph.name", "Row.names","ENTREZ","ACCNUM","SYMBOL","DESC"))
  write.csv(selected.reshape,file=paste0(working.dir,"/OUTPUT/selected.tableu.aug.plier.txt"), row.names = T)
}


#################################### ANNOTATION ######################################
all.data =1
if(all.data==1)
{
  # table.for.diff            <- core.data[core.data$ENTREZ!="NA",]
  table.for.diff            <- core.data
  table.for.diff$ACCNUM     <- NULL
  table.for.diff$DESC       <- NULL
  # gene.names                <- factor(paste(table.for.diff$SYMBOL,row.names(table.for.diff),"at",sep = "_"))
  # row.names(table.for.diff) <- gene.names
  table.for.diff$SYMBOL     <- NULL
  names(table.for.diff)     <- all.groups
  p.distribution            <- NULL
  t.distribution            <- NULL
  
  curr.diff             <- as.data.frame(table.for.diff[,1:96])
  g                     <- as.factor(c(rep("SalinePCB",12),rep("SalineLBS",12), rep("MagnesiumPCB",12), rep("MagnesiumLBS",12),
                             rep("BetaPCB",12), rep("BetaLBS",12), rep("MagBetaPCB",12),rep("MagBetaLBS",12)))
  design                <- model.matrix(~0+g)
  colnames(design)         <- levels(g)
  # colnames(curr.diff)      <- groups
  fit                   <- lmFit(curr.diff,design) 
  cont.matrix           <- makeContrasts(SalinePCB-SalineLBS,MagnesiumPCB-MagnesiumLBS,BetaPCB-BetaLBS,MagBetaPCB-MagBetaLBS, levels=design)
  fit2                  <- contrasts.fit(fit,cont.matrix)
  fit2                  <- eBayes(fit2)
  topTable(fit2,coef=1, adjust="bonferroni")
  graph.name            <- paste("All Groups Plier") 
  
  
  selected              <- p.adjust(fit2$p.value < 0.01, method = "bonferroni")
  row_distance          <- dist(as.matrix(curr.diff[selected==1,]), method = "manhattan")
  row_cluster           <- hclust(row_distance)
  col_distance          <- dist(t(as.matrix(curr.diff[selected==1,])), method = "manhattan")
  col_cluster           <- hclust(col_distance, method = "ward.D")
  
  png(file=paste0(working.dir,"/graph/allgenes-",graph.name,".png"),width = 900, height = 900, units = "px")
  heatmap.2(as.matrix(curr.diff), main = graph.name, key=TRUE, symkey=FALSE, density.info="none", trace="none", scale="none"
             ,ColSideColors = c(rep("green", 12),rep("blue", 12), rep("yellow",12),rep("red",12),rep("black",12),rep("brown",12),rep("white",12),rep=("darkorchid")),na.color = "black")
}
