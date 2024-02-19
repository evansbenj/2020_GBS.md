# AngSD for association tests

Path
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2022_allofraseri/mapped_to_XLv10_concatscaf
```

Example sbatch script:
```
#!/bin/sh
#SBATCH --job-name=angsd_allo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=256gb
#SBATCH --output=angsd_allo.%J.out
#SBATCH --error=angsd_allo.%J.err
#SBATCH --account=def-ben


module load StdEnv/2020 angsd/0.939

angsd -yBin bin_sex.ybin -doAsso 1 -GL 1 -out out_additive_F1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minInd 4 -ba
m bam.filelist -P 5 -doCounts 1 -setMinDepthInd 2 -setMaxDepthInd 100 -Pvalue 1
```
Filter to remove non-significant sites for plotting. 
Without the "-Pvalue 1' flag, the 6th column is a chisq value with df=1, so lets save only significant values (higher than 7):
```
zcat tempty.lrt0.gz | awk '$6 < 7 { next } { print }'> sig.only
zcat out_additive_F1.lrt0.gz | awk '$7 > 0.001 { next } { print }'> 2022_pyg_P_gt_0.001_only.txt
```

Plot whole genome and individual chromosomes:
```R
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/mapped_to_XL_v10_only")
library(tidyverse)
#install.packages("qqman")
library(qqman)
library(dplyr)



### BELOW works using lattice graphics instead of ggplot
# but there are issues with adding rectangles, etc
# but otherwise it looks pretty good
library(lattice)
library(plyr)
library(latticeExtra)
# manhattan plot function ----
# This script is from:
# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab=expression(-log[10](p)),
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(-log10(pvalue),thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}

# chromosome plot function ----

chromosome.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab=expression(-log[10](p)),
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fonte=18, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(-log10(pvalue),thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 #at=((posmax+posmin)/2+posshift),
                 labels=T, 
                 ticks=T, rot=0,
                 text.cex=1,
                 #tck = seq(0,200,50),
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    } else if(side=="left") {
      panel.axis(side=side, outside=T,
                 #at=((posmax+posmin)/2+posshift),
                 labels=T, 
                 ticks=T, rot=0,
                 text.cex=1,
                 #tck = seq(0,200,50),
                 check.overlap=F
      )
    } else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}

theme.novpadding <- list(
  layout.heights = list(
    top.padding = 0,
    main.key.padding = 0,
    key.axis.padding = 0,
    axis.xlab.padding = 0,
    xlab.key.padding = 0,
    key.sub.padding = 0,
    bottom.padding = 0
  ),
  layout.widths = list(
    left.padding = 0,
    key.ylab.padding = 0,
    ylab.axis.padding = 0,
    axis.key.padding = 0,
    right.padding = 0
  )
)

# muelleri ----
# read the data from muelleri
#my_df <- read.table("muelleri_additive_lt_0.1_only.txt", header = F)
my_df <- read.table(gzfile("muel_XL_v10_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
muel_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                       "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                       "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                       "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                                        my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                        "SIG"=list(col="red",label=F)),
                                        ylim= c(0,10),xlab = "",
                                        par.settings = theme.novpadding,
                                        panel=function(x, y, ...){
                                        panel.rect(xleft=695, ybottom=0,
                                        xright=733, ytop=10, alpha=0.3, col="light blue")
                                        panel.text(275,9,labels=expression(italic("X. muelleri")),fontsize=14)})
muel_wgs_plot


#pdf("./muelleri_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
jpeg("./muelleri_manhattan_angsd.jpg",w=14, h=3.0, units ="in", bg="transparent", res = 200)
  muel_wgs_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr4L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

#p2 <- xyplot(ts[12], type = "p", col = "red", cex = 10)

muel_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                    SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                                                ylim = c(0,8),
                                                xlab = " ",ylab = "",
                                                par.settings = theme.novpadding,
                                                panel=function(x, y, ...){
                                                panel.rect(xleft=111, ybottom=0,
                                                xright=147, ytop=8, alpha=0.3, col="light blue")
                                                panel.rect(xleft=36.4, ybottom=0,
                                                xright=36.6, ytop=8, col="black")
                                                panel.text(40,7,labels=expression(paste(italic("X. muelleri")," Chr4L")),fontsize=20)})



muel_sexchr_plot


#pdf("./muelleri_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
jpeg("./muelleri_manhattan_angsd_SexChr_only.jpg",w=7, h=3.0, units ="in", bg="transparent", res = 200)
  muel_sexchr_plot
dev.off()

# fischbergi ----
# read the data from fisch
my_df <- read.table(gzfile("fisc_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
#my_df[my_df$pvalue == 1,] <-0.9
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
fisc_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                                                ylim= c(0,17),xlab="",
                                                par.settings = theme.novpadding,
                                                panel=function(x, y, ...){
                                                panel.rect(xleft=468, ybottom=0,
                                                xright=532, ytop=17, alpha=0.3, col="light blue")
                                                panel.text(295,15,labels=expression(italic("X. fischbergi")),fontsize=14)})
fisc_wgs_plot

#pdf("./fisc_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
jpeg("./fisc_manhattan_angsd.jpg",w=14, h=3.0, units ="in", bg="transparent", res = 200)
  fisc_wgs_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr3L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


fisc_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)),
                                                  ylim = c(0,17),
                                                  xlab = "",
                                                  par.settings = theme.novpadding,
                                                  panel=function(x, y, ...){
                                                  panel.rect(xleft=41, ybottom=0,
                                                  xright=105, ytop=17, alpha=0.3, col="light blue")
                                                  panel.rect(xleft=19.1, ybottom=0,
                                                  xright=19.3, ytop=17, col="black")
                                                  panel.text(40,15,labels=expression(paste(italic("X. fischbergi")," Chr3L")),fontsize=20)})

fisc_sexchr_plot 


#pdf("./fisch_manhattan_angsd_Chr3L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
jpeg("./fisc_manhattan_angsd_SexChr_only.jpg",w=7, h=3.0, units ="in", bg="transparent", res = 200)
  fisc_sexchr_plot
dev.off()

# fraseri ----
# read the data from fisch
my_df <- read.table(gzfile("fraseri_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

#make annotation factor
#ann<-rep(1, length(my_df$pvalue))
#ann[with(my_df, pvalue<=0.0001)]<-2
#ann[with(my_df, pvalue<=0.00001)]<-3
#ann[with(my_df, chr=="Chr2L" & pos>=108800000 & pos<109400000)]<-4
#ann[with(my_df, chr=="Chr3L" & pos>=19600000 & pos<20800000)]<-5
#ann[with(my_df, chr=="Chr5S" & pos>=51600000 & pos<53000000)]<-6
#ann<-factor(ann, levels=1:6, labels=c("NS","MID","SIG","chr2L","chr3L","chr5S"))




# make the manhattan plot using the custom function above
fras_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                    ylim= c(0,10),xlab="",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                      panel.text(237,8,labels=expression(italic("X. fraseri")),fontsize=14)})


pdf("./fras_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  fras_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr3L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./fras_manhattan_angsd_Chr3L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  p_SL
dev.off()


# parafraseri ----
# read the data from fisch
my_df <- read.table(gzfile("parafras_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

# make the manhattan plot using the custom function above
para_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),                                                
                    ylim= c(0,10),xlab="",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                       panel.text(283,8,labels=expression(italic("X. parafraseri")),fontsize=14)})


pdf("./parafras_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
para_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr2S",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./parafras_manhattan_angsd_notreallySexChr2S_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# wittei ----
# read the data from wittei
my_df <- read.table(gzfile("wittei_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
witt_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                    ylim= c(0,12),xlab = "",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                      panel.text(220,10,labels=expression(italic("X. wittei")),fontsize=14)})


pdf("./witt_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr3L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./witt_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# longipes ----
# read the data from longipes
my_df <- read.table(gzfile("long_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

# make the manhattan plot using the custom function above
long_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                    ylim= c(0,14),xlab = "",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                      panel.text(261,9,labels=expression(italic("X. longipes")),fontsize=14)})


pdf("./long_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
long_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr9_10S",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./long_manhattan_angsd_notreallySexChr9_10S_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
    p_SL
dev.off()

# lenduensis ----
# read the data from lenduensis
my_df <- read.table(gzfile("lend_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

# make the manhattan plot using the custom function above
lend_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                                                ylim= c(0,17),xlab="",
                                                par.settings = theme.novpadding,
                                                panel=function(x, y, ...){
                                                panel.rect(xleft=440, ybottom=0,
                                                xright=441.2, ytop=17, alpha=0.3, col="light blue")
                                                panel.text(302,15,labels=expression(italic("X. lenduensis")),fontsize=14)})
                    


#pdf("./lend_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
jpeg("./lend_manhattan_angsd.jpg",w=14, h=3.0, units ="in", bg="transparent", res = 200)
  lend_wgs_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr3L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
#ann[with(SexChr, chr=="Chr3L" & pos>=15073758 & pos<16244640)]<-4
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))



lend_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)),
                                                    ylim = c(0,12),
                                                    xlab = "",ylab = "",
                                                    par.settings = theme.novpadding,
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=15, ybottom=0,
                                                    xright=16.2, ytop=12, alpha=0.3, col="light blue")
                                                    panel.rect(xleft=19.1, ybottom=0,
                                                    xright=19.3, ytop=12, col="black")
                                                    panel.text(43,10,labels=expression(paste(italic("X. lenduensis")," Chr3L")),fontsize=20)})


lend_sexchr_plot 

#pdf("./lend_manhattan_angsd_Chr3L__only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
jpeg("./lend_manhattan_angsd_SexChr_only.jpg",w=7, h=3.0, units ="in", bg="transparent", res = 200)
  lend_sexchr_plot
dev.off()



# largeni ----
# read the data from largeni
my_df <- read.table(gzfile("larg_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))


# make the manhattan plot using the custom function above
larg_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                 "SIG"=list(col="red",label=F),
                                                 "neardmw"=list(col="green",label=F),
                                                 "dmw"=list(col="blue",label=F)),                   
                    ylim= c(0,10),xlab = "",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                      panel.text(240,8,labels=expression(italic("X. largeni")),fontsize=14)})


pdf("./larg_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
larg_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr2L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann[with(SexChr, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(SexChr, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                   "SIG"=list(col="red",label=F),
                                                   "neardmw"=list(col="green",label=F),
                                                   "dmw"=list(col="blue",label=F)))


pdf("./larg_manhattan_angsd_notreallySexChr2L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# itombwensis ----
# read the data from itombwensis
my_df <- read.table(gzfile("itom_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))


# make the manhattan plot using the custom function above
itom_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F),
                                                "neardmw"=list(col="green",label=F),
                                                "dmw"=list(col="blue",label=F)),
                    ylim= c(0,12),xlab = "",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                      panel.text(297,10,labels=expression(italic("X. itombwensis")),fontsize=14)})


pdf("./itom_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
 itom_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr2L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann[with(SexChr, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(SexChr, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F),
                                                    "neardmw"=list(col="green",label=F),
                                                    "dmw"=list(col="blue",label=F)))


pdf("./itom_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# cliv_all ----
# read the data from cliv_all
my_df <- read.table(gzfile("cliv_all_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))


# make the manhattan plot using the custom function above
cliv_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F),
                                                "neardmw"=list(col="green",label=F),
                                                "dmw"=list(col="blue",label=F)),
                    ylim= c(0,15),xlab = "",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                      panel.text(211,13,labels=expression(italic("X. clivii")),fontsize=14)})



pdf("./cliv_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
cliv_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann[with(SexChr, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(SexChr, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))

p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F),
                                                    "neardmw"=list(col="green",label=F),
                                                    "dmw"=list(col="blue",label=F)))


pdf("./cliv_all_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  p_SL
dev.off()


# cliv_onlyEritrea ----
# read the data from cliv_onlyEritrea
my_df <- read.table(gzfile("clivonlyEritrea_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))

# make the manhattan plot using the custom function above
cliv_erit_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F),
                                                "neardmw"=list(col="green",label=F),
                                                "dmw"=list(col="blue",label=F)),
                    ylim= c(0,10),xlab = "",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                      panel.text(300,8,labels=expression(paste("Eritrea ",italic("X. clivii"))),fontsize=14)})


pdf("./clivonlyEritrea_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  cliv_erit_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr5L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./clivonlyEritrea_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# laevis lab ----
# read the data from laevis
my_df <- read.table(gzfile("laevis_out_Xlaev_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))



laev_lab_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                "SIG"=list(col="red",label=F),
                                                                "neardmw"=list(col="green",label=F),
                                                                "dmw"=list(col="blue",label=F)),
                                                    ylim= c(0,11),xlab="",
                                                    par.settings = theme.novpadding,
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=405, ybottom=0,
                                                    xright=425, ytop=11, alpha=0.3, col="light blue")
                                                    panel.text(280,10,labels=expression(paste("lab ",italic("X. laevis"))),fontsize=14)})




#pdf("./laev_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
jpeg("./laev_family_manhattan_angsd.jpg",w=14, h=3.0, units ="in", bg="transparent", res = 200)
  laev_lab_wgs_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr2L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
#ann[with(SexChr, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
#ann[with(SexChr, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))#,"neardmw","dmw"))



laev_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)),
                                                    #"neardmw"=list(col="green",label=F),
                                                    #"dmw"=list(col="blue",label=F)),
                                                    ylim = c(0,11),
                                                    xlab = "",
                                                    par.settings = theme.novpadding,
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=177, ybottom=0,
                                                    xright=190, ytop=11, alpha=0.3, col="light blue")
                                                    panel.rect(xleft=70.1, ybottom=0,
                                                    xright=70.3, ytop=11, col="black")
                                                    panel.text(40,9,labels=expression(paste(italic("X. laevis")," Chr2L")),fontsize=20)})

laev_sexchr_plot 

#p_SL
#pdf("./laev_family_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
jpeg("./laev_fam_manhattan_angsd_SexChr_only.jpg",w=7, h=3.0, units ="in", bg="transparent", res = 200)
  laev_sexchr_plot
dev.off()

# wild laevis ----
# read the data from wild laevis
my_df <- read.table(gzfile("XLwild_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
my_df[my_df$pvalue == 0,]$pvalue <-0.0000000000000032196
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))




# make the manhattan plot using the custom function above
wild_xl_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                              my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                          "SIG"=list(col="red",label=F),
                                                          "neardmw"=list(col="green",label=F),
                                                          "dmw"=list(col="blue",label=F)),
                              ylim= c(0,17),xlab="",
                              par.settings = theme.novpadding,
                              panel=function(x, y, ...){
                                panel.rect(xleft=405, ybottom=0,
                                           xright=425, ytop=17, alpha=0.3, col="light blue")
                                panel.text(300,15,labels=expression(paste("wild ",italic("X. laevis"))),fontsize=14)})

wild_xl_wgs
pdf("./wildlaev_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
wild_xl_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr2L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
#ann[with(SexChr, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
#ann[with(SexChr, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))#,"neardmw","dmw"))




p_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)),
                                                    #"neardmw"=list(col="green",label=F),
                                                    #"dmw"=list(col="blue",label=F))),
                                                    ylim = c(0,17),
                                                    xlab = "Chr2L",
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=182, ybottom=0,
                                                    xright=183, ytop=17, alpha=0.3, col="light blue")
                                                    panel.rect(xleft=70.1, ybottom=0,
                                                    xright=70.3, ytop=17, alpha=0.3, col="black")
                                                    panel.text(40,15,labels=expression(paste("wild ",italic("X. laevis"))),fontsize=20)})

p_SL

pdf("./wildlaev_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  p_SL
dev.off()


# vict ----
# read the data from vict
my_df <- read.table(gzfile("vict_out_Xlaev_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))


# make the manhattan plot using the custom function above
vict_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F),
                                                "neardmw"=list(col="green",label=F),
                                                "dmw"=list(col="blue",label=F)),
                                                ylim= c(0,10),xlab = "",
                                                par.settings = theme.novpadding,
                                                panel=function(x, y, ...){
                                                panel.text(275,9,labels=expression(italic("X. victorianus")),fontsize=14)})

pdf("./vict_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  vict_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr2L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann[with(SexChr, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(SexChr, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))



p_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                    "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                    "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                    "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                        SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                     "SIG"=list(col="red",label=F),
                                                     "neardmw"=list(col="green",label=F),
                                                     "dmw"=list(col="blue",label=F)),
                                                      #"neardmw"=list(col="green",label=F),
                                                      #"dmw"=list(col="blue",label=F))),
                                                      ylim = c(0,17),
                                                      xlab = "Chr2L",
                                                      panel=function(x, y, ...){
                                                      panel.rect(xleft=182, ybottom=0,
                                                      xright=183, ytop=17, alpha=0.3, col="light blue")
                                                      panel.rect(xleft=70.1, ybottom=0,
                                                      xright=70.3, ytop=17, alpha=0.3, col="black")
                                                      panel.text(40,15,labels=expression(italic("X. victorianus")),fontsize=20)})

p_SL


pdf("./vict_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  p_SL
dev.off()

# vict uncertain----
# read the data from vict including 3 uncertain sexes with preliminary designations
my_df <- read.table(gzfile("vict_plusuncertain_out_Xlaev_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=0.0001)]<-2
ann[with(my_df, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))

pdf("./vict_uncertain_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr3L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./vict_uncertain_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# gilli ----
# read the data from gilli
my_df <- read.table(gzfile("gilli_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))


# make the manhattan plot using the custom function above
gill_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                      "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                      "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                      "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                           my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                       "SIG"=list(col="red",label=F),
                                                       "neardmw"=list(col="green",label=F),
                                                       "dmw"=list(col="blue",label=F)),
                           ylim= c(0,10),xlab = "",
                           par.settings = theme.novpadding,
                           panel=function(x, y, ...){
                             panel.text(208,8,labels=expression(italic("X. gilli")),fontsize=14)})



#gilli_wgs

pdf("./gill_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
gill_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr2L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
#ann[with(SexChr, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
#ann[with(SexChr, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))#,"neardmw","dmw"))


p_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)),
                                                    #"neardmw"=list(col="green",label=F),
                                                    #"dmw"=list(col="blue",label=F))),
                                                    ylim = c(0,5),
                                                    xlab = "Chr2L",
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=182, ybottom=0,
                                                    xright=183, ytop=5, alpha=0.3, col="light blue")
                                                    panel.rect(xleft=70.1, ybottom=0,
                                                    xright=70.3, ytop=5, alpha=0.3, col="black")
                                                    panel.text(40,4,labels=expression(italic("X. gilli")),fontsize=20)})


p_SL

pdf("./gill_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  p_SL
dev.off()


# tropicalis - Ghana ----
# read the data from trop_all - this is only from Ghana not from Mitros
my_df <- read.table(gzfile("out_Xtrop_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
# rename chromosomes
my_df[my_df$chr == "Chr1:1-217471166",]$chr <- "Chr1"
my_df[my_df$chr == "Chr2:1-181034961",]$chr <- "Chr2"
my_df[my_df$chr == "Chr3:1-153873357",]$chr <- "Chr3"
my_df[my_df$chr == "Chr4:1-153961319",]$chr <- "Chr4"
my_df[my_df$chr == "Chr5:1-164033575",]$chr <- "Chr5"
my_df[my_df$chr == "Chr6:1-154486312",]$chr <- "Chr6"
my_df[my_df$chr == "Chr7:1-133565930",]$chr <- "Chr7"
my_df[my_df$chr == "Chr8:1-147241510",]$chr <- "Chr8"
my_df[my_df$chr == "Chr9:1-91218944",]$chr <- "Chr9"
my_df[my_df$chr == "Chr10:1-52432566",]$chr <- "Chr10"
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scafs",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
trop_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                               "Chr5","Chr6","Chr7","Chr8","Chr9",
                                               "Chr10")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                                                ylim= c(0,17),xlab = "",
                                                par.settings = theme.novpadding,
                                                panel=function(x, y, ...){
                                                panel.rect(xleft=1024, ybottom=0,
                                                xright=1040, ytop=17, alpha=0.3, col="light blue")
                                                panel.text(190,15,labels=expression(paste(italic("X. tropicalis")," Ghana")),
                                                fontsize=14)})

pdf("./trop_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
    trop_wgs_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

trop_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                   "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                   "Chr10")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)),
                                                    ylim = c(0,17),
                                                    xlab = "",ylab = "",
                                                    par.settings = theme.novpadding,
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=6, ybottom=0,
                                                    xright=11, ytop=17, alpha=0.3, col="light blue")
                                                    panel.rect(xleft=60.3, ybottom=0,
                                                    xright=60.5, ytop=17, col="black")
                                                    panel.text(45,14,labels=expression(paste(italic("X. tropicalis")," Ghana Chr7")),fontsize=20)})

trop_sexchr_plot 


pdf("./trop_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  trop_sexchr_plot
dev.off()

# tropicalis - GW ----
# read the data from trop_GW
my_df <- read.table(gzfile("GW_out_Xtrop_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-15
# rename chromosomes
my_df[my_df$chr == "Chr1:1-217471166",]$chr <- "Chr1"
my_df[my_df$chr == "Chr2:1-181034961",]$chr <- "Chr2"
my_df[my_df$chr == "Chr3:1-153873357",]$chr <- "Chr3"
my_df[my_df$chr == "Chr4:1-153961319",]$chr <- "Chr4"
my_df[my_df$chr == "Chr5:1-164033575",]$chr <- "Chr5"
my_df[my_df$chr == "Chr6:1-154486312",]$chr <- "Chr6"
my_df[my_df$chr == "Chr7:1-133565930",]$chr <- "Chr7"
my_df[my_df$chr == "Chr8:1-147241510",]$chr <- "Chr8"
my_df[my_df$chr == "Chr9:1-91218944",]$chr <- "Chr9"
my_df[my_df$chr == "Chr10:1-52432566",]$chr <- "Chr10"
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scafs",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                               "Chr5","Chr6","Chr7","Chr8","Chr9",
                                               "Chr10")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))

pdf("./trop_GW_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                   "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                   "Chr10")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./trop_GW_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# tropicalis - GE ----
# read the data from trop_GE
my_df <- read.table(gzfile("GE_out_Xtrop_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-15
# rename chromosomes
my_df[my_df$chr == "Chr1:1-217471166",]$chr <- "Chr1"
my_df[my_df$chr == "Chr2:1-181034961",]$chr <- "Chr2"
my_df[my_df$chr == "Chr3:1-153873357",]$chr <- "Chr3"
my_df[my_df$chr == "Chr4:1-153961319",]$chr <- "Chr4"
my_df[my_df$chr == "Chr5:1-164033575",]$chr <- "Chr5"
my_df[my_df$chr == "Chr6:1-154486312",]$chr <- "Chr6"
my_df[my_df$chr == "Chr7:1-133565930",]$chr <- "Chr7"
my_df[my_df$chr == "Chr8:1-147241510",]$chr <- "Chr8"
my_df[my_df$chr == "Chr9:1-91218944",]$chr <- "Chr9"
my_df[my_df$chr == "Chr10:1-52432566",]$chr <- "Chr10"
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scafs",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                               "Chr5","Chr6","Chr7","Chr8","Chr9",
                                               "Chr10")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))

pdf("./trop_GE_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                   "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                   "Chr10")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./trop_GE_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# tropicalis - wildonly ----
# read the data from trop_wildonly
my_df <- read.table(gzfile("Xtrop_wild_caught_only_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-15
# rename chromosomes
my_df[my_df$chr == "Chr1:1-217471166",]$chr <- "Chr1"
my_df[my_df$chr == "Chr2:1-181034961",]$chr <- "Chr2"
my_df[my_df$chr == "Chr3:1-153873357",]$chr <- "Chr3"
my_df[my_df$chr == "Chr4:1-153961319",]$chr <- "Chr4"
my_df[my_df$chr == "Chr5:1-164033575",]$chr <- "Chr5"
my_df[my_df$chr == "Chr6:1-154486312",]$chr <- "Chr6"
my_df[my_df$chr == "Chr7:1-133565930",]$chr <- "Chr7"
my_df[my_df$chr == "Chr8:1-147241510",]$chr <- "Chr8"
my_df[my_df$chr == "Chr9:1-91218944",]$chr <- "Chr9"
my_df[my_df$chr == "Chr10:1-52432566",]$chr <- "Chr10"
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scafs",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
trop_wild_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                               "Chr5","Chr6","Chr7","Chr8","Chr9",
                                               "Chr10")), my_df$pos, 
                                my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),
                                ylim= c(0,17),xlab = "",
                                par.settings = theme.novpadding,
                                panel=function(x, y, ...){
                                  panel.rect(xleft=1024, ybottom=0,
                                             xright=1040, ytop=17, alpha=0.3, col="light blue")
                                  panel.text(168,15,labels=expression(paste("wild ",italic("X. tropicalis"))),
                                             fontsize=14)})

pdf("./trop_wildonly_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  trop_wild_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                   "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                   "Chr10")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./trop_wildonly_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# mellotropicalis ----
# read the data from mello
my_df <- read.table(gzfile("out_Xmello_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-15
# rename chromosomes
my_df[my_df$chr == "Chr1:1-217471166",]$chr <- "Chr1"
my_df[my_df$chr == "Chr2:1-181034961",]$chr <- "Chr2"
my_df[my_df$chr == "Chr3:1-153873357",]$chr <- "Chr3"
my_df[my_df$chr == "Chr4:1-153961319",]$chr <- "Chr4"
my_df[my_df$chr == "Chr5:1-164033575",]$chr <- "Chr5"
my_df[my_df$chr == "Chr6:1-154486312",]$chr <- "Chr6"
my_df[my_df$chr == "Chr7:1-133565930",]$chr <- "Chr7"
my_df[my_df$chr == "Chr8:1-147241510",]$chr <- "Chr8"
my_df[my_df$chr == "Chr9:1-91218944",]$chr <- "Chr9"
my_df[my_df$chr == "Chr10:1-52432566",]$chr <- "Chr10"
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scafs",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))



# make the manhattan plot using the custom function above
mell_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                               "Chr5","Chr6","Chr7","Chr8","Chr9",
                                               "Chr10")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                    ylim= c(0,9),xlab = "",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                      panel.rect(xleft=1024, ybottom=0,
                                 xright=1040, ytop=17, alpha=0.3, col="light blue")
                      panel.text(178,7,labels=expression(italic("X. mellotropicalis")),fontsize=14)})



pdf("./mello_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
mell_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

p_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                   "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                   "Chr10")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)),
                                                    ylim = c(0,6),
                                                    xlab = "Chr7",
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=0.25, ybottom=0,
                                                    xright=9.5, ytop=6, alpha=0.3, col="light blue")})
p_SL



pdf("./mello_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  p_SL
dev.off()


# boum ----
# read the data from boum
my_df <- read.table(gzfile("boum_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=0.0001)]<-2
ann[with(my_df, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
boum_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                                        ylim= c(0,10),xlab = "",
                                        par.settings = theme.novpadding,
                                       panel=function(x, y, ...){
                                       panel.text(311,8,labels=expression(italic("X. boumbaensis")),fontsize=14)})


pdf("./boum_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
    boum_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr2L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./boum_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# borealis lab family ----
# read the data from borealis
my_df <- read.table(gzfile("bor_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
bore_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)),
                                                ylim= c(0,17),xlab="",
                                                par.settings = theme.novpadding,
                                                panel=function(x, y, ...){
                                                panel.rect(xleft=1210, ybottom=0,
                                                xright=1290, ytop=17, alpha=0.3, col="light blue")
                                                panel.text(315,15,labels=expression(paste("lab ",italic("X. borealis"))),fontsize=14)})
bore_wgs_plot
#pdf("./bore_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
jpeg("./bore_manhattan_angsd.jpg",w=14, h=3.0, units ="in", bg="transparent", res = 200)
  bore_wgs_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr8L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


bore_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)),
                                                    ylim = c(0,17),
                                                    xlab = "",
                                                    par.settings = theme.novpadding,
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=0, ybottom=0,
                                                    xright=71.5, ytop=17, alpha=0.3, col="light blue")
                                                    panel.rect(xleft=22, ybottom=0,
                                                    xright=22.3, ytop=17, col="black")
                                                    panel.text(43,15,labels=expression(paste(italic("X. borealis")," east Chr8L",sep=" ")),
                                                               fontsize=20)})

bore_sexchr_plot 


pdf("./bore_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  bore_sexchr_plot
dev.off()


# wild borealis west ----
# read the data from borealis
my_df <- read.table(gzfile("XB_wild_westonly_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
wild_bor_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                                my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),
                                ylim= c(0,12),xlab="",
                                par.settings = theme.novpadding,
                                panel=function(x, y, ...){
                                  panel.rect(xleft=1210, ybottom=0,
                                             xright=1290, ytop=17, alpha=0.3, col="light blue")
                                  panel.text(336,10,labels=expression(paste(italic("X. borealis")," west",sep=" ")),fontsize=14)})
wild_bor_wgs_plot
#pdf("./bore_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
jpeg("./wild_bore_west_manhattan_angsd.jpg",w=14, h=3.0, units ="in", bg="transparent", res = 200)
  wild_bor_wgs_plot
dev.off()

# wild borealis east ----
# read the data from borealis
my_df <- read.table(gzfile("XB_wild_eastonly_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

# make the manhattan plot using the custom function above
wild_bor_east_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                "SIG"=list(col="red",label=F)),
                                    ylim= c(0,10),xlab="",
                                    par.settings = theme.novpadding,
                                    panel=function(x, y, ...){
                                      panel.rect(xleft=1210, ybottom=0,
                                                 xright=1290, ytop=17, alpha=0.3, col="light blue")
                                      panel.text(336,8,labels=expression(paste(italic("X. borealis")," east",sep=" ")),fontsize=14)})
wild_bor_east_wgs_plot
#pdf("./bore_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
jpeg("./wild_bore_east_manhattan_angsd.jpg",w=14, h=3.0, units ="in", bg="transparent", res = 200)
  wild_bor_east_wgs_plot
dev.off()



pdf("./bore_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
wild_bor_sexchr_plot
dev.off()


# allo_family1 ----
# read the data from allo_fam1
my_df <- read.table(gzfile("allo_family1_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=0.0001)]<-2
ann[with(my_df, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=0.0001)]<-2
ann[with(my_df, pvalue<=0.00001)]<-3
ann[with(my_df, chr=="Chr7L" & pos>=18481359 & pos<19973768)]<-4
ann<-factor(ann, levels=1:4, labels=c("NS","MID","SIG","SL"))


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F),
                                                "SL"=list(col="green",label=F)))

pdf("./allo_fam1_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann[with(SexChr, chr=="Chr7L" & pos>=18481359 & pos<19973768)]<-4
ann<-factor(ann, levels=1:4, labels=c("NS","MID","SIG","SL"))

p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                   "SIG"=list(col="red",label=F),
                                                   "SL"=list(col="green",label=F)))


pdf("./allo_fam1_manhattan_angsd_Chr7L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# allo_family2 ----
# read the data from allo_fam2
my_df <- read.table(gzfile("allo_family2_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=0.005)]<-2
ann[with(my_df, pvalue<=0.0001)]<-3
ann[with(my_df, pvalue<=0.00001)]<-4
ann<-factor(ann, levels=1:4, labels=c("NS","MID","SIG","VSIG"))


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="lightblue",label=F),
                                                "SIG"=list(col="orange",label=F), 
                                                "VSIG"=list(col="red",label=F)))

pdf("./allo_fam2_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann[with(SexChr, chr=="Chr7L" & pos>=18481359 & pos<19973768)]<-4
ann<-factor(ann, levels=1:4, labels=c("NS","MID","SIG","SL"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F),
                                                    "SL"=list(col="green",label=F)))


pdf("./allo_fam2_manhattan_angsd_Chr7_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# allo_all ----
# read the data from allo_all
my_df <- read.table(gzfile("all_allo_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
my_df[my_df$pvalue == 0,]$pvalue <-2.3315e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

#ann<-rep(1, length(my_df$pvalue))
#ann[with(my_df, pvalue<=0.0001)]<-2
#ann[with(my_df, pvalue<=0.00001)]<-3
#ann[with(my_df, chr=="Chr7L" & pos>=18481359 & pos<19973768)]<-4
#ann<-factor(ann, levels=1:4, labels=c("NS","MID","SIG","SL"))

# make the manhattan plot using the custom function above
allo_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                              "SIG"=list(col="red",label=F)),
                                              ylim= c(0,17),xlab="",
                                              par.settings = theme.novpadding,
                                              panel=function(x, y, ...){
                                              panel.rect(xleft=1088, ybottom=0,
                                              xright=1105, ytop=17, alpha=0.3, col="light blue")
                                              panel.text(295,15,labels=expression(italic("X. allofraseri")),fontsize=14)})


allo_wgs_plot

pdf("./allo_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  allo_wgs_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

#make annotation factor
#ann<-rep(1, length(SexChr$pvalue))
#ann[with(SexChr, pvalue<=0.0001)]<-2
#ann[with(SexChr, pvalue<=0.00001)]<-3
#ann[with(SexChr, chr=="Chr7L" & pos>=18481359 & pos<19973768)]<-4
#ann<-factor(ann, levels=1:4, labels=c("NS","MID","SIG","SL"))

allo_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                   "SIG"=list(col="red",label=F)),
                                                    ylim= c(0,17),
                                                    xlab = "",
                                                    par.settings = theme.novpadding,
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=8, ybottom=0,
                                                    xright=25, ytop=17, alpha=0.3, col="light blue")
                                                    panel.rect(xleft=58.2, ybottom=0,
                                                    xright=58.5, ytop=17, col="black")
                                                    panel.text(35,15,labels=expression(paste(italic("X. allofraseri")," Chr7L")),fontsize=20)})

allo_sexchr_plot 

pdf("./allo_all_manhattan_angsd_Chr7L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  allo_sexchr_plot
dev.off()



# pygm_all ----
# read the data from pygm_all
my_df <- read.table(gzfile("pygm_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))



# make the manhattan plot using the custom function above
pygm_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                                my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F),
                                                            "neardmw"=list(col="green",label=F),
                                                            "dmw"=list(col="blue",label=F)),
                                ylim= c(0,16),xlab="",
                                                par.settings = theme.novpadding,
                                                panel=function(x, y, ...){
                                                panel.rect(xleft=1333, ybottom=0,
                                                xright=1355, ytop=16, alpha=0.3, col="light blue")
                                                panel.text(300,14,labels=expression(paste(italic("X. pygmaeus"))),fontsize=14)})

pygm_wgs_plot

pdf("./pygm_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  pygm_wgs_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr8L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))



pygm_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)),
                                                    ylim= c(0,17),
                                                    xlab = "",ylab = "",
                                                    par.settings = theme.novpadding,
                                                    panel=function(x, y, ...){
                                                    panel.rect(xleft=117, ybottom=0,
                                                    xright=135.4, ytop=17, alpha=0.3, col="light blue")
                                                    panel.rect(xleft=22, ybottom=0,
                                                    xright=22.3, ytop=17, col="black")
                                                    panel.text(35,14,labels=expression(paste(italic("X. pygmaeus")," Chr8L")),fontsize=20)})


pygm_sexchr_plot 


pdf("./pygm_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  pygm_sexchr_plot
dev.off()


# pygm_lab ----
# read the data from pygm_lab
my_df <- read.table(gzfile("pygm_only_labreared_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))



# make the manhattan plot using the custom function above
pygm_lab_wgs_plot <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                                my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F),
                                                            "neardmw"=list(col="green",label=F),
                                                            "dmw"=list(col="blue",label=F)),
                                ylim= c(0,16),xlab="",
                                par.settings = theme.novpadding,
                                panel=function(x, y, ...){
                                  panel.rect(xleft=1333, ybottom=0,
                                             xright=1355, ytop=16, alpha=0.3, col="light blue")
                                  panel.text(340,14,labels=expression(paste("lab ",italic("X. pygmaeus"))),fontsize=14)})

pygm_lab_wgs_plot

pdf("./pygm_lab_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  pygm_lab_wgs_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr8L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))



pygm_labonly_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                                "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                                "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                                "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                                    SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                 "SIG"=list(col="red",label=F)),
                                    ylim= c(0,17),
                                    xlab = "",ylab = "",
                                    par.settings = theme.novpadding,
                                    panel=function(x, y, ...){
                                      panel.rect(xleft=117, ybottom=0,
                                                 xright=135.4, ytop=17, alpha=0.3, col="light blue")
                                      panel.rect(xleft=22, ybottom=0,
                                                 xright=22.3, ytop=17, col="black")
                                      panel.text(35,14,labels=expression(paste("lab ",italic("X. pygmaeus")," Chr8L")),fontsize=20)})


pygm_labonly_sexchr_plot 


pdf("./pygm_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
pygm_labonly_sexchr_plot
dev.off()


# pygm_wildonly ----
# read the data from pygm_wildonly
my_df <- read.table(gzfile("pygm_only_wild_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann[with(my_df, chr=="Chr2L" & pos>=180693575 & pos<184720555)]<-4
ann[with(my_df, chr=="Chr2L" & pos>=182693575 & pos<182720555)]<-5
ann<-factor(ann, levels=1:5, labels=c("NS","MID","SIG","neardmw","dmw"))


# make the manhattan plot using the custom function above
wild_pygm_wgs <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                                my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F),
                                                            "neardmw"=list(col="green",label=F),
                                                            "dmw"=list(col="blue",label=F)),
                                ylim= c(0,10),xlab="",
                    par.settings = theme.novpadding,
                    panel=function(x, y, ...){
                      panel.rect(xleft=1333, ybottom=0,
                                 xright=1355, ytop=16, alpha=0.3, col="light blue")
                      panel.text(350,8,labels=expression(paste("wild ",italic("X. pygmaeus"))),fontsize=14)})

pdf("./pygm_wildonly_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
wild_pygm_wgs
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr8L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


wild_pygm_sexchr_plot <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                                "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                                "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                                "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                                    SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                 "SIG"=list(col="red",label=F)),
                                                                ylim= c(0,10),
                                                                xlab = "Chr8L",ylab = "",
                                                                par.settings = theme.novpadding,
                                                                panel=function(x, y, ...){
                                                                panel.rect(xleft=117, ybottom=0,
                                                                xright=135.4, ytop=10, alpha=0.3, col="light blue")
                                                                panel.rect(xleft=22, ybottom=0,
                                                                xright=22.3, ytop=10, col="black")
                                                                panel.text(35,8,labels=expression(paste("wild ",italic("X. pygmaeus"))),fontsize=20)})


pdf("./pygm_wildonly_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
wild_pygm_sexchr_plot
dev.off()






# MAT PAT only ---- 
# for Xl
# Chr1L	233740090
# Chr1S	202412970
# Chr2L	191000146
# Chr2S	169306100
# Chr3L	161426101
# Chr3S	131962816
# Chr4L	155250554
# Chr4S	132731174
# Chr5L	171415384
# Chr5S	143394103
# Chr6L	164223595
# Chr6S	137316286
# Chr7L	139837618
# Chr7S	113060389
# Chr8L	135449133
# Chr8S	103977862
# Chr9_10L	137811819
# Chr9_10S	117266291

# allo_fam0_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("allo_fam0_mat_angsdout.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))



# make the manhattan plot using the custom function above
allo_fam0_mat_all_plot <- chromosome.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))




pdf("./allo_fam0_mat_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  allo_fam0_mat_all_plot
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

allomat_fam0_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                                     "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                                     "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                                     "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                                   SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                      "SIG"=list(col="red",label=F)), 
                                                                      ylim = c(0,5),
                                                                      xlab = "",
                                                                      par.settings = theme.novpadding,
                                                                      panel=function(x, y, ...){
                                                                      panel.rect(xleft=8, ybottom=0,
                                                                      xright=25, ytop=5, alpha=0.3, col="light blue")
                                                                      panel.rect(xleft=58.2, ybottom=0,
                                                                      xright=58.5, ytop=5, col="black")})


pdf("./allo_fam0_mat_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  allomat_fam0_SL
dev.off()

# allo_fam0_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("allo_fam0_pat_angsdout.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
allopat_fam0_SL_all <- chromosome.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))


pdf("./allo_fam0_pat_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  allopat_fam0_SL_all
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))




allopat_fam0_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                    "SIG"=list(col="red",label=F)), 
                                                                    ylim = c(0,5),
                                                                    xlab = "",
                                                                   par.settings = theme.novpadding,
                                                                   panel=function(x, y, ...){
                                                                   panel.rect(xleft=8, ybottom=0,
                                                                   xright=25, ytop=5, alpha=0.3, col="light blue")
                                                                  panel.rect(xleft=58.2, ybottom=0,
                                                                  xright=58.5, ytop=5, col="black")})

pdf("./allo_fam0_pat_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  allopat_fam0_SL
dev.off()

library(gridExtra)
pdf("./allo_fam0_matpat_angsd_Chr7L_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(allomat_fam0_SL, allopat_fam0_SL, ncol=2, 
               nrow =1, 
               top = "X. allofraseri family0")
              #top = textGrob(expression(paste(italic("X. allofraseri")," family0",sep=" ")),gp=gpar(fontsize=20,font=3)))
dev.off()

# allo_fam1_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("allo_family1_mat_angsdout", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))



# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))

pdf("./allo_fam1_mat_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

allomat_fam1_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                         "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                         "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                         "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                             SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                          "SIG"=list(col="red",label=F)), 
                                                          ylim = c(0,15),
                                                          xlab = "",
                                                          par.settings = theme.novpadding,
                                                          panel=function(x, y, ...){
                                                          panel.rect(xleft=8, ybottom=0,
                                                          xright=25, ytop=17, alpha=0.3, col="light blue")
                                                          panel.rect(xleft=58.2, ybottom=0,
                                                          xright=58.5, ytop=17, col="black")})

allomat_fam1_SL 


pdf("./allo_fam1_mat_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  allomat_fam1_SL
dev.off()

# allo_fam1_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("allo_family1_pat_angsdout", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))

pdf("./allo_fam1_pat_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


allopat_fam1_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                              "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                              "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                              "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                                  SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                               "SIG"=list(col="red",label=F)),
                                                                ylab="", ylim = c(0,15),
                                                                xlab = "",                              
                                                                par.settings = theme.novpadding,
                                                                panel=function(x, y, ...){
                                                                panel.rect(xleft=8, ybottom=0,
                                                                xright=25, ytop=17, alpha=0.3, col="light blue")
                                                                panel.rect(xleft=58.2, ybottom=0,
                                                                xright=58.5, ytop=17, col="black")})

allopat_fam1_SL 


pdf("./allo_fam1_pat_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  allopat_fam1_SL
dev.off()

library(gridExtra)


pdf("./allo_fam1_matpat_angsd_Chr7L_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(allomat_fam1_SL, allopat_fam1_SL, ncol=2, nrow =1)
dev.off()


# allo_fam2_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("allo_fam2_mat_angsdout.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))



# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))

pdf("./allo_fam2_mat_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


allomat_fam2_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                               "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                                   SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                "SIG"=list(col="red",label=F)), 
                                                                ylim = c(0,5),
                                                                xlab = "",
                                                                par.settings = theme.novpadding,
                                                                panel=function(x, y, ...){
                                                                panel.rect(xleft=8, ybottom=0,
                                                                xright=25, ytop=5, alpha=0.3, col="light blue")
                                                                panel.rect(xleft=58.2, ybottom=0,
                                                                xright=58.5, ytop=5, col="black")})

pdf("./allo_fam2_mat_all_manhattan_angsd_Chr7L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  allomat_fam2_SL
dev.off()

# allo_fam2_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("allo_fam2_pat_angsdout.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scaffolds",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=a["0.1%"])]<-2
ann[with(my_df, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))

pdf("./allo_fam2_pat_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


allopat_fam2_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                               "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                                   SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                "SIG"=list(col="red",label=F)), 
                                                                ylim = c(0,5),
                                                                xlab = "",
                                                                par.settings = theme.novpadding,
                                                                panel=function(x, y, ...){
                                                                panel.rect(xleft=8, ybottom=0,
                                                                xright=25, ytop=5, alpha=0.3, col="light blue")
                                                                panel.rect(xleft=58.2, ybottom=0,
                                                                xright=58.5, ytop=5, col="black")})


pdf("./allo_fam2_pat_all_manhattan_angsd_Chr7L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  allopat_fam2_SL
dev.off()

library(gridExtra)


pdf("./allo_fam2_matpat_angsd_Chr7L_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(allomat_fam2_SL, allopat_fam2_SL, ncol=2, nrow =1)
dev.off()

# muel_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("muel_XL_v10_out_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr4L
SexChr <- read.table("muel_Chr4L_family1_mat_angsdout", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


muel_mat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                               "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                                   SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                "SIG"=list(col="red",label=F)),xlab="", ylim = c(0,8),
                                                                par.settings = theme.novpadding,
                                                                panel=function(x, y, ...){
                                                                panel.rect(xleft=111, ybottom=0,
                                                                xright=147, ytop=8, alpha=0.3, col="light blue")
                                                                panel.rect(xleft=36.4, ybottom=0,
                                                                xright=36.7, ytop=8, col="black")
                                                                panel.text(38,6,labels=expression(paste(italic("X. muelleri")," Chr4L")),fontsize=14)})

muel_mat_SL 

pdf("./muel_fam1_mat_all_manhattan_angsd_Chr4L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
muel_mat_SL
dev.off()

# muel_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("muel_XL_v10_out_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr4L
SexChr <- read.table("muel_Chr4L_family1_pat_angsdout", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


muel_pat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),xlab="", 
                                                            ylab="", ylim = c(0,8),
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=111, ybottom=0,
                                                            xright=147, ytop=8, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=36.4, ybottom=0,
                                                            xright=36.7, ytop=8, col="black")})

muel_pat_SL 



pdf("./muel_fam1_pat_all_manhattan_angsd_Chr4L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
muel_pat_SL
dev.off()

library(gridExtra)


pdf("./muel_fam1_matpat_angsd_Chr4L_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(muel_mat_SL, muel_pat_SL, ncol=2, nrow =1)
dev.off()





# fisc_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("fisc_out_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr3L
SexChr <- read.table("fisc_Chr3L_family1_mat_angsdout", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

fisc_mat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),xlab="", ylim = c(0,8),
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=41, ybottom=0,
                                                            xright=105, ytop=8, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=19.1, ybottom=0,
                                                            xright=19.4, ytop=8, col="black")
                                                            panel.text(42,6.5,labels=expression(paste(italic("X. fischbergi")," Chr3L")),fontsize=14)})

fisc_mat_SL

pdf("./fisc_fam1_mat_all_manhattan_angsd_Chr4L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
fisc_mat_SL
dev.off()

# fisc_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("fisc_out_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr3L
SexChr <- read.table("fisc_Chr3L_family1_pat_angsdout", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


fisc_pat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),xlab="",
                                                            ylab="",ylim = c(0,8),
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=41, ybottom=0,
                                                            xright=105, ytop=8, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=19.1, ybottom=0,
                                                            xright=19.4, ytop=8, col="black")})

fisc_pat_SL

pdf("./fisc_fam1_pat_all_manhattan_angsd_Chr3L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  fisc_pat_SL
dev.off()

library(gridExtra)
pdf("./fisc_fam1_matpat_angsd_Chr3L_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
grid.arrange(fisc_mat_SL, fisc_pat_SL, ncol=2, nrow =1)
dev.off()




# bore_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("bor_out_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16

# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr8L
SexChr <- read.table("bore_Chr8L_mat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

bore_mat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),xlab="", 
                                                            xlim = c(0,135.449133),
                                                            ylim = c(0,11),
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=0, ybottom=0,
                                                            xright=71.5, ytop=12, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=22, ybottom=0,
                                                            xright=22.3, ytop=12, col="black")
                                                            panel.text(33,9,labels=expression(paste(italic("X. borealis")," Chr8L")),fontsize=14)})

bore_mat_SL

pdf("./bore_pat_mat_all_manhattan_angsd_Chr8L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  bore_mat_SL
dev.off()

# bore_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("bor_out_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr8L
SexChr <- read.table("bore_Chr8L_pat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


bore_pat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),xlab="", 
                                                            xlim = c(0,135.449133),
                                                            ylab="", ylim = c(0,11),
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=0, ybottom=0,
                                                            xright=71.5, ytop=12, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=22, ybottom=0,
                                                            xright=22.3, ytop=12, col="black")})

bore_pat_SL

pdf("./bore_pat_all_manhattan_angsd_Chr8L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  bore_pat_SL
dev.off()

library(gridExtra)
pdf("./bore_pat_SL_matpat_angsd_Chr8L_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(bore_mat_SL, bore_pat_SL, ncol=2, nrow =1)
dev.off()




# pygm_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("pygm_only_labreared_out_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr8L
SexChr <- read.table("pyg_Chr8L_mat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

pygm_mat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),xlab="", ylim = c(0,17),
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=117, ybottom=0,
                                                            xright=135.4, ytop=17, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=22, ybottom=0,
                                                            xright=22.3, ytop=17, col="black")
                                                            panel.text(35,13,labels=expression(paste(italic("X. pygmaeus")," Chr8L")),fontsize=14)})

pygm_mat_SL

pdf("./pygm_pat_mat_all_manhattan_angsd_Chr8L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
pygm_mat_SL
dev.off()

# pygm_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("pygm_only_labreared_out_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr8L
SexChr <- read.table("pyg_Chr8L_pat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


pygm_pat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),xlab="", 
                                                            ylab="", ylim = c(0,17),
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=117, ybottom=0,
                                                            xright=135.4, ytop=17, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=22, ybottom=0,
                                                            xright=22.3, ytop=17, col="black")})

pygm_pat_SL 

pdf("./pygm_pat_all_manhattan_angsd_Chr8L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  pygm_pat_SL
dev.off()

library(gridExtra)
pdf("./pygm_pat_SL_matpat_angsd_Chr8L_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(pygm_mat_SL, pygm_pat_SL, ncol=2, nrow =1)
dev.off()


# laev_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("laevis_out_Xlaev_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr8L
SexChr <- read.table("XL_Chr2L_mat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

laev_mat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),xlab="", ylim = c(0,12),
                                                              par.settings = theme.novpadding,
                                                              panel=function(x, y, ...){
                                                              panel.rect(xleft=177, ybottom=0,
                                                              xright=190, ytop=12, alpha=0.3, col="light blue")
                                                              panel.rect(xleft=70.1, ybottom=0,
                                                              xright=70.3, ytop=12, col="black")
                                                              panel.text(40,9,labels=expression(paste(italic("X. laevis")," Chr2L")),fontsize=14)})

laev_mat_SL

pdf("./laev_pat_mat_all_manhattan_angsd_Chr2L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  laev_mat_SL
dev.off()

# laev_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("laevis_out_Xlaev_additive_F1.lrt0.gz", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

# now open only mat sites for Chr8L
SexChr <- read.table("XL_Chr2L_pat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["0.1%"])]<-2
ann[with(SexChr, pvalue<=a["0.05%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


laev_pat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                           "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                           "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                           "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)),xlab="", 
                                                            ylab="", ylim = c(0,12),
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=177, ybottom=0,
                                                            xright=190, ytop=12, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=70.1, ybottom=0,
                                                            xright=70.3, ytop=12, col="black")})

laev_pat_SL

pdf("./laev_pat_all_manhattan_angsd_Chr2L_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  laev_pat_SL
dev.off()

library(gridExtra)
pdf("./laev_pat_SL_matpat_angsd_Chr8L_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
grid.arrange(laev_mat_SL, laev_pat_SL, ncol=2, nrow =1)
dev.off()



# trop_Mitros_C659_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("Mitros_trop_C659_mat_angsdout.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.01, 0.005))

# now open only C659 mat sites for Chr7
SexChr <- read.table("Mitros_trop_C659_mat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["1%"])]<-2
ann[with(SexChr, pvalue<=a["0.5%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


C659_mat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                              "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                              "Chr10")), SexChr$pos, 
                                  SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                               "SIG"=list(col="red",label=F)), 
                                                              ylim = c(0,17),
                                                              xlab = "",
                                                              par.settings = theme.novpadding,
                                                              panel=function(x, y, ...){
                                                              panel.rect(xleft=6, ybottom=0,
                                                              xright=11, ytop=17, alpha=0.3, col="light blue")
                                                              panel.rect(xleft=60.3, ybottom=0,
                                                              xright=60.6, ytop=17, col="black")})

pdf("./C659_mat_manhattan_angsd_Chr7_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  C659_mat_SL
dev.off()

# trop_Mitros_C659_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("Mitros_trop_C659_pat_angsdout.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.01, 0.005))

# now open only C659 mat sites for Chr7
SexChr <- read.table("Mitros_trop_C659_pat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["1%"])]<-2
ann[with(SexChr, pvalue<=a["0.5%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

C659_pat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                           "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                           "Chr10")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)), 
                                                            ylim = c(0,17),
                                                            xlab = "",
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=6, ybottom=0,
                                                            xright=11, ytop=17, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=60.3, ybottom=0,
                                                            xright=60.6, ytop=17, col="black")})

pdf("./C659_pat_manhattan_angsd_Chr7_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  C659_pat_SL
dev.off()

library(gridExtra)
pdf("./C659_matpat_angsd_Chr7_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(C659_mat_SL, C659_pat_SL, ncol=2, nrow =1)
dev.off()




# trop_Mitros_C660_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("Mitros_trop_C660_mat_angsdout.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.01, 0.005))

# now open only C660 mat sites for Chr7
SexChr <- read.table("Mitros_trop_C660_mat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["1%"])]<-2
ann[with(SexChr, pvalue<=a["0.5%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

C660_mat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                           "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                           "Chr10")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)), 
                                                            ylim = c(0,12),
                                                            xlab = "",
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=6, ybottom=0,
                                                            xright=11, ytop=12, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=60.3, ybottom=0,
                                                            xright=60.6, ytop=12, col="black")})



pdf("./C660_mat_manhattan_angsd_Chr7_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  C660_mat_SL
dev.off()

# trop_Mitros_C660_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("Mitros_trop_C660_pat_angsdout.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.01, 0.005))

# now open only C660 mat sites for Chr7
SexChr <- read.table("Mitros_trop_C660_pat_angsdout.txt", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["1%"])]<-2
ann[with(SexChr, pvalue<=a["0.5%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


C660_pat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                           "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                           "Chr10")), SexChr$pos, 
                               SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                            "SIG"=list(col="red",label=F)), 
                                                            ylim = c(0,12),
                                                            xlab = "",
                                                            par.settings = theme.novpadding,
                                                            panel=function(x, y, ...){
                                                            panel.rect(xleft=6, ybottom=0,
                                                            xright=11, ytop=12, alpha=0.3, col="light blue")
                                                            panel.rect(xleft=60.3, ybottom=0,
                                                            xright=60.6, ytop=12, col="black")})


pdf("./C660_pat_manhattan_angsd_Chr7_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  C660_pat_SL
dev.off()

library(gridExtra)
pdf("./C660_matpat_angsd_Chr7_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(C660_mat_SL, C660_pat_SL, ncol=2, nrow =1)
dev.off()





# trop_GE_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("GE_trop__mat_angsdout", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.01, 0.005))

# now open only C660 mat sites for Chr7
SexChr <- read.table("GE_trop__mat_angsdout", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["1%"])]<-2
ann[with(SexChr, pvalue<=a["0.5%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


GE_trop_mat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                              "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                              "Chr10")), SexChr$pos, 
                                  SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                               "SIG"=list(col="red",label=F)), 
                                                                ylim = c(0,4),
                                                                xlab = "",
                                                                par.settings = theme.novpadding,
                                                                panel=function(x, y, ...){
                                                                panel.rect(xleft=6, ybottom=0,
                                                                xright=11, ytop=4, alpha=0.3, col="light blue")
                                                                panel.rect(xleft=60.3, ybottom=0,
                                                                xright=60.6, ytop=4, col="black")})


pdf("./GE_trop_mat_manhattan_angsd_Chr7_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  GE_trop_mat_SL
dev.off()

# GE_trop_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("GE_trop__pat_angsdout", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.01, 0.005))

# now open only C660 mat sites for Chr7
SexChr <- read.table("GE_trop__pat_angsdout", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["1%"])]<-2
ann[with(SexChr, pvalue<=a["0.5%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


GE_trop_pat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                              "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                              "Chr10")), SexChr$pos, 
                                  SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                               "SIG"=list(col="red",label=F)), 
                                                                ylim = c(0,4),
                                                                xlab = "",
                                                                par.settings = theme.novpadding,
                                                                panel=function(x, y, ...){
                                                                panel.rect(xleft=6, ybottom=0,
                                                                xright=11, ytop=4, alpha=0.3, col="light blue")
                                                                panel.rect(xleft=60.3, ybottom=0,
                                                                xright=60.6, ytop=4, col="black")})


pdf("./GE_trop_pat_manhattan_angsd_Chr7_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  GE_trop_pat_SL
dev.off()

library(gridExtra)
pdf("./GE_trop_matpat_angsd_Chr7_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(GE_trop_mat_SL, GE_trop_pat_SL, ncol=2, nrow =1)
dev.off()







# trop_GW_mat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("GW_trop__mat_angsdout", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.01, 0.005))

# now open only GW mat sites for Chr7
SexChr <- read.table("GW_trop__mat_angsdout", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["1%"])]<-2
ann[with(SexChr, pvalue<=a["0.5%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

GW_trop_mat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                              "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                              "Chr10")), SexChr$pos, 
                                   SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                                "SIG"=list(col="red",label=F)), 
                                                                ylim = c(0,5),
                                                                xlab = "",
                                                                par.settings = theme.novpadding,
                                                                panel=function(x, y, ...){
                                                                panel.rect(xleft=6, ybottom=0,
                                                                xright=11, ytop=5, alpha=0.3, col="light blue")
                                                                panel.rect(xleft=60.3, ybottom=0,
                                                                xright=60.6, ytop=5, col="black")})

pdf("./GW_trop_mat_manhattan_angsd_Chr7_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  GW_trop_mat_SL
dev.off()

# GW_trop_pat_only ----
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/gets_matpat')
my_df <- read.table("GW_trop__pat_angsdout", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
# get percentiles from whole genome including all sites
a <- quantile(my_df$pvalue, probs = c(0.01, 0.005))

# now open only GW pat sites for Chr7
SexChr <- read.table("GW_trop__pat_angsdout", header = T)
colnames(SexChr) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

#View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
# Get rid of the scaffold data

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=a["1%"])]<-2
ann[with(SexChr, pvalue<=a["0.5%"])]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


GW_trop_pat_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                              "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                              "Chr10")), SexChr$pos, 
                                  SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                               "SIG"=list(col="red",label=F)), 
                                                               ylim = c(0,5),
                                                               xlab = "",
                                                               par.settings = theme.novpadding,
                                                               panel=function(x, y, ...){
                                                               panel.rect(xleft=6, ybottom=0,
                                                               xright=11, ytop=5, alpha=0.3, col="light blue")
                                                               panel.rect(xleft=60.3, ybottom=0,
                                                               xright=60.6, ytop=5, col="black")})


pdf("./GW_trop_pat_manhattan_angsd_Chr7_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  GW_trop_pat_SL
dev.off()

library(gridExtra)
pdf("./GW_trop_matpat_angsd_Chr7_only.pdf",w=7, h=2.0, version="1.4", bg="transparent")
  grid.arrange(GW_trop_mat_SL, GW_trop_pat_SL, ncol=2, nrow =1)
dev.off()




# tropicalis - WGS no tads ----
# read the data from trop_all
setwd('/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/mapped_to_XL_v10_only')
my_df <- read.table("2023_XT_notads_WGS_P_gt_0.001_only.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero

#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scafs",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=0.0001)]<-2
ann[with(my_df, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                               "Chr5","Chr6","Chr7","Chr8","Chr9",
                                               "Chr10")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)), ylim = c(3,6));p

pdf("./XT_notads_WGS_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr5",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

p_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                   "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                   "Chr10")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)), ylim = c(3,6));p_SL


pdf("./trop_notads_WGS_angsd_Chr5_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  p_SL
dev.off()


# tropicalis - WGS with tads ----
# read the data from trop_all
setwd('/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/mapped_to_XL_v10_only')
my_df <- read.table("XT_withtads_P_gt_0.001_only.txt", header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
#View(my_df) # check if there are any pvalues equal to zero
my_df[my_df$chr == "Chr1:1-217471166",]$chr <- "Chr1"
my_df[my_df$chr == "Chr2:1-181034961",]$chr <- "Chr2"
my_df[my_df$chr == "Chr3:1-153873357",]$chr <- "Chr3"
my_df[my_df$chr == "Chr4:1-153961319",]$chr <- "Chr4"
my_df[my_df$chr == "Chr5:1-164033575",]$chr <- "Chr5"
my_df[my_df$chr == "Chr6:1-154486312",]$chr <- "Chr6"
my_df[my_df$chr == "Chr7:1-133565930",]$chr <- "Chr7"
my_df[my_df$chr == "Chr8:1-147241510",]$chr <- "Chr8"
my_df[my_df$chr == "Chr9:1-91218944",]$chr <- "Chr9"
my_df[my_df$chr == "Chr10:1-52432566",]$chr <- "Chr10"
#my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
# Get rid of the scaffold data
my_df <- my_df[my_df$chr != "Scafs",]

# get percentiles
a <- quantile(my_df$pvalue, probs = c(0.001, 0.0005))

#make annotation factor
ann<-rep(1, length(my_df$pvalue))
ann[with(my_df, pvalue<=0.0001)]<-2
ann[with(my_df, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                               "Chr5","Chr6","Chr7","Chr8","Chr9",
                                               "Chr10")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)), ylim = c(3,8));p

pdf("./XT_with_tads_WGS_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))

p_SL <- chromosome.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                    "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                    "Chr10")), SexChr$pos, 
                        SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                     "SIG"=list(col="red",label=F)), 
                        xlim = c(0,15), ylim = c(3,8));p_SL


pdf("./trop_with_tads_WGS_angsd_Chr7_0_15Mb_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# multiple Manhattan plots ----
library(gridExtra)

# Fig 2 ----
jpeg("./Fig2_multiple_sexchronly_manhattan_angsd.jpg",w=12, h=8.0, units ="in", bg="transparent", res = 200)
      grid.arrange(laev_sexchr_plot,pygm_sexchr_plot,
               allo_sexchr_plot,lend_sexchr_plot,
               fisc_sexchr_plot,muel_sexchr_plot, 
               bore_sexchr_plot,trop_sexchr_plot, ncol=2, nrow =4)
dev.off()

# SI Fig1a ---- 
jpeg("./SI_Fig1a_multiple_sexchronly_manhattan_angsd.jpg",w=12, h=13.0, units ="in", bg="transparent", res = 200)
      grid.arrange(laev_lab_wgs_plot,wild_xl_wgs,
                   pygm_lab_wgs_plot,wild_pygm_wgs,
                   bore_wgs_plot,wild_bor_east_wgs_plot,
                   wild_bor_wgs_plot,ncol=1, nrow =8)
dev.off()

# SI Fig1b ---- 
jpeg("./SI_Fig1b_multiple_wgs_noSL_manhattan_angsd.jpg",w=12, h=13.0, units ="in", bg="transparent", res = 200)
      grid.arrange(cliv_wgs,cliv_erit_wgs,
             vict_wgs,gill_wgs,
             larg_wgs,allo_wgs_plot,
             lend_wgs_plot,fisc_wgs_plot,ncol=1, nrow =8)
dev.off()

# SI Fig1c ---- 
jpeg("./SI_Fig1c_multiple_wgs_noSL_manhattan_angsd.jpg",w=12, h=13.0, units ="in", bg="transparent", res = 200)
    grid.arrange(fras_wgs,para_wgs,
             muel_wgs_plot,itom_wgs,
             boum_wgs,witt_wgs,
             long_wgs,ncol=1, nrow=8)
dev.off()

# SI Fig1d ---- 
jpeg("./SI_Fig1d_multiple_wgs_noSL_manhattan_angsd.jpg",w=12, h=13.0, units ="in", bg="transparent", res = 200)
    grid.arrange(trop_wgs_plot,trop_wild_wgs,
             mell_wgs,ncol=1, nrow =8)
dev.off()


library(gridExtra)
library(ggpubr)
jpeg("./allo_matpat_family0_and_2_angsd.jpg",w=10, h=3.0, units ="in", bg="transparent", res = 200)
  grid.arrange(arrangeGrob(allomat_fam0_SL,allomat_fam2_SL, top = "Maternal"),
             arrangeGrob(allopat_fam0_SL,allopat_fam2_SL, top = "Paternal"), ncol = 2)
dev.off()

# Fig3 ---- 
title1=text_grob("Maternal                                                              Paternal", size = 16)
jpeg("./Fig3_allo_matpat_family1_angsd.jpg",w=10, h=3, units ="in", bg="transparent", res = 200)
      grid.arrange(arrangeGrob(allomat_fam1_SL),
      arrangeGrob(allopat_fam1_SL), ncol = 2, top=title1)
dev.off()

# SI Fig3 ---- 
title1=text_grob("Maternal                                                              Paternal", size = 16)
jpeg("./SI_Fig3_female_heterogamy_matpat_5_species_angsd.jpg",w=10, h=8.0, units ="in", bg="transparent", res = 200)
  grid.arrange(laev_mat_SL, laev_pat_SL,
             pygm_mat_SL,pygm_pat_SL,
             bore_mat_SL,bore_pat_SL,
             fisc_mat_SL,fisc_pat_SL,
             muel_mat_SL, muel_pat_SL, ncol = 2, nrow=5, top=title1)
dev.off()


# SI Fig4 ---- 
title1=text_grob("Maternal                                                              Paternal", size = 16)
jpeg("./SI_Fig4_allo_matpat_2families_angsd.jpg",w=12, h=4.0, units ="in", bg="transparent", res = 200)
grid.arrange(allomat_fam0_SL, allopat_fam0_SL, 
             allomat_fam2_SL, allopat_fam2_SL, ncol=2, nrow =2, top=title1)
dev.off()

# SI Fig5 ---- 
title1=text_grob("Maternal                                                              Paternal", size = 16)
jpeg("./SI_Fig5_trop_matpat_4_families_angsd.jpg",w=10, h=8.0, units ="in", bg="transparent", res = 200)
grid.arrange(GE_trop_mat_SL, GE_trop_pat_SL,
             GW_trop_mat_SL,GW_trop_pat_SL,
             C659_mat_SL,C659_pat_SL,
             C660_mat_SL,C660_pat_SL, ncol = 2, nrow =4, top=title1)
dev.off()

```
