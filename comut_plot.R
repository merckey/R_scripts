# example for co-mut plot 
# per http://stackoverflow.com/questions/25300822/r-plot-to-show-mutation
# with heavy modification and hacking
# read ggplot2 from Hadley Wickham for a lot of details 
# good explanation of geom_raster:
# http://rgm3.lab.nig.ac.jp/RGM/R_rdfile?f=ggplot2/man/geom_raster.Rd&d=R_CC
# geom_raster() is much faster than geom_tile()

library(ggplot2)

# data frame for mutation profile of gene, subj, mut 
# mut will be represented by NA, 0, 1, etc.  - can change to nonsynominous, etc
# geom_raster will be plotted for "mut" with 0 or 1 value but not NA
# I changed this to saying missense, nonsense, transversion to actually simulate real data
# lists 10 genes
# with 50 subjects

dat <- expand.grid(gene=1:10, subj=1:50)
dat$mut <- as.factor(sample(c(rep("Missense",300),rep("nonsense",200),rep("transversion",10)),500))
dat$mut[sample(500,300)] <- NA

# can't use geom_raster as can't set size of boxes
mut <- ggplot(dat, aes(x=subj, y=gene, height=0.8, width=0.8))
mut + geom_tile(aes(fill=mut)) + 
  #scale_fill_manual(values = c("red","orange","black"), na.value="Grey90") +
  scale_fill_brewer(palette = "Set1", na.value="Grey90") +
  scale_x_continuous(breaks=1:50) + 
  xlab("Subject") +
  scale_y_continuous(breaks=1:10,labels=c("D0","D1","D2","D3","D4","D5",
                                          "D6","Brca2","p53","Kras")) +
  guides(fill=F) +
  theme(
    axis.ticks.x=element_blank(), 
    axis.ticks.y=element_blank(),
    axis.text.x=element_blank(), 
    axis.text.y=element_text(colour="Black"), 
    axis.title.x=element_text(face="bold"), 
    axis.title.y=element_blank(),
    panel.grid.major.x=element_blank(), 
    panel.grid.major.y=element_blank(),
    panel.grid.minor.x=element_blank(), 
    panel.grid.minor.y=element_blank(), 
    panel.background=element_blank()
    )

# to set scale to -log10 for FDR 
# scale_x_continuous(trans = "-log10")
# or 
# p <- qplot(-log10(FDR), sample, data=comut_data)
# p + xlab("-log10 q-values")

#save plot data to disk 
save(mut2, file = "practice_comutplot2.rdata")
# save plot to disk
ggsave("practice_comutplot2.png",width=5, height = 3)




# alternative
ggplot(dat, aes(x=subj, y=gene, fill=mut)) + 
  geom_raster() +
  scale_fill_manual(values = c("red","orange","black"), na.value="#FFFFFF")

sessionInfo()
