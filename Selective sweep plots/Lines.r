# User: Shuo Liu
# Name: liushuo028@163.com
# Date: 2020.5-2020.11
# Purpose of the scripts: plot the selective sweep signals


library(ggplot2)
ggplot(dat, aes(pos_start, Value, color = group, group = group)) +     ## use one selective sweep result table (example.csv)
  geom_line(aes(linetype=group))+									
  scale_linetype_manual(values=c("twodash", "dotted","solid"))+		   ## specify the types of lines
  xlim(13100000, 13510000)+											   
  ylim(0, 0.012)+													   ## set x,y axis scaling
  xlab("Chromosome 4") + 
  ylab(expression(paste("Diversity (", italic("Pi"),")")))+
  theme_classic()+theme(legend.position='none')+
  scale_color_manual(values=c("black", "black", "black"))