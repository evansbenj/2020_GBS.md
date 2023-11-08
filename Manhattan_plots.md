# Makin g manhattan plots with angsd output
```
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

# This script is from:
# https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R
manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab=expression(-log[10](p-value)),
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

# muelleri ----
# read the data from muelleri
#my_df <- read.table("muelleri_additive_lt_0.1_only.txt", header = F)
my_df <- read.table(gzfile("muel_XL_v10_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")

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

pdf("./muelleri_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
  p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr4L",]

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


pdf("./muelleri_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  p_SL
dev.off()

# fischbergi ----
# read the data from fisch
my_df <- read.table(gzfile("fisch_XLv10_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
#my_df[my_df$pvalue == 1,] <-0.9
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

pdf("./fisch_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
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


pdf("./fisch_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
  p_SL
dev.off()

# fraseri ----
# read the data from fisch
my_df <- read.table(gzfile("fraseri_out_additive_F1.lrt0.gz"), header = T)
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

pdf("./fras_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
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


pdf("./fras_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# parafraseri ----
# read the data from fisch
my_df <- read.table(gzfile("parafras_out_additive_F1.lrt0.gz"), header = T)
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

pdf("./parafras_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr8L",]

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


pdf("./parafras_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# wittei ----
# read the data from wittei
my_df <- read.table(gzfile("wittei_out_additive_F1.lrt0.gz"), header = T)
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

pdf("./witt_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
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


pdf("./witt_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# longipes ----
# read the data from longipes
my_df <- read.table(gzfile("long_out_additive_F1.lrt0.gz"), header = T)
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

pdf("./long_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
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


pdf("./long_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# lenduensis ----
# read the data from lenduensis
my_df <- read.table(gzfile("lend_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
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

pdf("./lend_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
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


pdf("./lend_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()



# largeni ----
# read the data from lenduensis
my_df <- read.table(gzfile("larg_out_additive_F1.lrt0.gz"), header = T)
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

pdf("./larg_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
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


pdf("./larg_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# itombwensis ----
# read the data from itombwensis
my_df <- read.table(gzfile("itom_out_additive_F1.lrt0.gz"), header = T)
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

pdf("./itom_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
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


pdf("./itom_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# cliv_all ----
# read the data from cliv_all
my_df <- read.table(gzfile("cliv_all_out_additive_F1.lrt0.gz"), header = T)
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

pdf("./cliv_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

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


pdf("./cliv_all_manhattan_angsd_notreallySexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# cliv_onlyEritrea ----
# read the data from cliv_onlyEritrea
my_df <- read.table(gzfile("clivonlyEritrea_out_additive_F1.lrt0.gz"), header = T)
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

pdf("./clivonlyEritrea_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
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


# laevis ----
# read the data from laevis
my_df <- read.table(gzfile("laevis_out_Xlaev_additive_F1.lrt0.gz"), header = T)
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

pdf("./laev_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
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


pdf("./laev_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# vict ----
# read the data from vict
my_df <- read.table(gzfile("vict_out_Xlaev_additive_F1.lrt0.gz"), header = T)
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

pdf("./vict_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr9_10S",]

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

# tropicalis - all ----
# read the data from trop_all
my_df <- read.table(gzfile("out_Xtrop_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
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
                                                "SIG"=list(col="red",label=F)))

pdf("./trop_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                   "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                   "Chr10")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./trop_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
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
                                                "SIG"=list(col="red",label=F)))

pdf("./trop_GW_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
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
View(my_df) # check if there are any pvalues equal to zero
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
                                                "SIG"=list(col="red",label=F)))

pdf("./trop_GE_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
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
View(my_df) # check if there are any pvalues equal to zero
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
                                                "SIG"=list(col="red",label=F)))

pdf("./trop_wildonly_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
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
View(my_df) # check if there are any pvalues equal to zero
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
                                                "SIG"=list(col="red",label=F)))

pdf("./mello_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7",]

#make annotation factor
ann<-rep(1, length(SexChr$pvalue))
ann[with(SexChr, pvalue<=0.0001)]<-2
ann[with(SexChr, pvalue<=0.00001)]<-3
ann<-factor(ann, levels=1:3, labels=c("NS","MID","SIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1","Chr2","Chr3","Chr4",
                                                   "Chr5","Chr6","Chr7","Chr8","Chr9",
                                                   "Chr10")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                    "SIG"=list(col="red",label=F)))


pdf("./mello_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# boum ----
# read the data from boum
my_df <- read.table(gzfile("boum_out_additive_F1.lrt0.gz"), header = T)
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

pdf("./boum_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
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

# borealis ----
# read the data from borealis
my_df <- read.table(gzfile("bor_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
my_df[my_df$pvalue == 0,]$pvalue <-1.1102e-16
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

pdf("./bore_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr8L",]

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


pdf("./bore_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
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


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))

pdf("./allo_fam1_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

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


pdf("./allo_fam1_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
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
ann[with(SexChr, pvalue<=0.001)]<-2
ann[with(SexChr, pvalue<=0.0001)]<-3
ann[with(SexChr, pvalue<=0.00001)]<-4
ann<-factor(ann, levels=1:4, labels=c("NS","MID","SIG","VSIG"))


p_SL <- manhattan.plot(factor(SexChr$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                                   "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                                   "Chr7S","Chr8S","Chr9_10S")), SexChr$pos, 
                       SexChr$pvalue, annotate=list(ann, "MID"=list(col="lightblue",label=F), 
                                                    "SIG"=list(col="orange",label=F),
                                                    "VSIG"=list(col="red",label=F)))


pdf("./allo_fam2_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# allo_all ----
# read the data from allo_all
my_df <- read.table(gzfile("all_allo_out_additive_F1.lrt0.gz"), header = T)
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


# make the manhattan plot using the custom function above
p <- manhattan.plot(factor(my_df$chr, levels=c("Chr1L","Chr2L","Chr3L","Chr4L",
                                               "Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                               "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S",
                                               "Chr7S","Chr8S","Chr9_10S")), my_df$pos, 
                    my_df$pvalue, annotate=list(ann, "MID"=list(col="orange",label=F), 
                                                "SIG"=list(col="red",label=F)))

pdf("./allo_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr7L",]

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


pdf("./allo_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# pygm_all ----
# read the data from pygm_all
my_df <- read.table(gzfile("pygm_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
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

pdf("./pygm_all_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr8L",]

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


pdf("./pygm_all_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

# pygm_wildonly ----
# read the data from pygm_wildonly
my_df <- read.table(gzfile("pygm_only_wild_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
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

pdf("./pygm_wildonly_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr8L",]

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


pdf("./pygm_wildonly_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()


# pygm_labonly ----
# read the data from pygm_labonly
my_df <- read.table(gzfile("pygm_only_labreared_out_additive_F1.lrt0.gz"), header = T)
colnames(my_df) <- c("chr","pos","MAJ","MIN","FREQ","LRT","pvalue")
View(my_df) # check if there are any pvalues equal to zero
#my_df[my_df$pvalue == 0,]$pvalue <-2.2204e-15
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

pdf("./pygm_labonly_manhattan_angsd.pdf",w=14, h=3.0, version="1.4", bg="transparent")
p
dev.off()

# Highlight Chromosome of interest
SexChr <- my_df[my_df$chr == "Chr8L",]

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


pdf("./pygm_labonly_manhattan_angsd_SexChr_only.pdf",w=7, h=3.0, version="1.4", bg="transparent")
p_SL
dev.off()

```
