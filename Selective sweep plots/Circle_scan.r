# User: Shuo Liu
# Name: liushuo028@163.com
# Date: 2020.5-2020.11
# Purpose of the scripts: plot the selective sweep signals in a circle

library(circlize)

dat <-read.csv("circos_data.csv",header = TRUE) ## selective sweep result table
ref<-read.table("Genome_len.chr",header = TRUE) ## length of chromosomes

## create an initial/basic circle
## set genome breaks, facing angle, fonts and so on
circos.clear()
col_text <- "grey20"
circos.par("track.height"=0.8,gap.degree=5,start.degree =86,clock.wise = T,
           cell.padding=c(0,0,0,0))
circos.initialize(factors=ref$Genome,										### see 
                  xlim=matrix(c(rep(0,8),ref$Length),ncol=2))
circos.track(ylim=c(0,1),panel.fun=function(x,y) {
  Genome=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim),mean(ylim),Genome,cex=0.5,col=col_text,
              facing="bending.inside",niceFacing=TRUE)
},bg.col="lightyellow2",bg.border=F,track.height=0.06)
brk <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5)*10^7
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h="top",major.at=brk,labels=round(brk/10^7,1),labels.cex=0.4,
              col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
},bg.border=F)

## an example: C1a_Likelihood.csv, and set threshold = 25.623 (dotted line)
circos.track(factors=dat$Chr,x=dat$pos_start,y=dat$C1a_Likelihood,
             panel.fun=function(x,y) {
               circos.lines(x,y,type="h", col="dimgrey",lwd=0.6)
               circos.segments(x0=0,x1=max(ref$Length),y0=0,y1=0,lwd=0.6,lty="11",col="dimgrey")
               circos.segments(x0=0,x1=max(ref$Length),y0=25.623,y1=25.623,lwd=0.6,lty="11",col="dimgrey")
             },ylim=c(0,87),track.height=0.08,bg.border=F)
circos.yaxis(side = c("right"),labels.cex=0.3,lwd=0.1,tick.length=0.1,labels.col=col_text,col="#FFFFFF")