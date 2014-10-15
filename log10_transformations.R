emt <- c("Epithelial",'Mesenchymal')
library(ggplot2)

setwd("~/Desktop/RNAseq_Nicole_Ecad/DESeq2_analysis_batch2/")
# did calculation in excel file 2014-10-9-cdh1-rpm.xls 
names <- c("PD2204",'PD2329','PD2342','PD2412','PD2523','PD345','PD798','PD9210')
rpm_ratio <- scan()
rpm_ratio
df <- data.frame(names, rpm_ratio)
cdh1_tab <- df

(logratio <- ggplot(cdh1_tab, aes(conds, log10(rpm_ratio))) + 
  geom_boxplot(outlier.color=NA) +
  geom_dotplot(binaxis='y', binwidth=0.025, stackdir='center',fill=NA))


(logratio <- logratio + ylab("log10 (RPM Ecad- / RPM Ecad+)"))
(logratio <- logratio + ggtitle("Cdh1 expression\np=2.3e-3"))
(logratio <- logratio + theme(
  plot.title = element_text(size = 20),
  axis.text.y=element_text(color='Black', face='bold',size=10), 
  axis.title.y=element_text(color='Black', face='bold',size=12), 
  axis.text.x=element_text(colour="Black",face="bold",size=10),
  axis.title.x=element_blank(),
  panel.grid.major.x=element_blank(), 
  panel.grid.major.y=element_blank(),
  panel.grid.minor.x=element_blank(), 
  panel.grid.minor.y=element_blank(), 
  panel.background=element_rect(color='black',fill=NA)
))

ggsave(sparc_ratio, width=3, height=5, file='2014-10-10-Sparc-rpm.png')
save(sparc_ratio, file="2014-10-10-ggplot2_log10targetlevels.RData")
save.image("2014-10-10-batch2.RData")

x<- assay(ddsMF['Sparc'])
write.csv(x,"~/Desktop/sparc.csv")
# added RPM ratio values for sparc from excel file sparc.csv
cdh1_tab$qpcr <- scan()

sparc_ratio <- ggplot(cdh1_tab, aes(conds, log10(sparc_rpmratio))) + geom_boxplot(outlier.color=NA)

(sparc_ratio <- sparc_ratio + ggtitle("Sparc expression\np=3.62e-2") +
  geom_dotplot(binaxis='y', binwidth=0.025, stackdir='center',fill=NA) +
   ylab("log10 (RPM Ecad- / RPM Ecad+)"+
  theme(
  plot.title = element_text(size = 20),
  axis.text.y=element_text(color='Black', face='bold',size=10), 
  axis.title.y=element_text(color='Black', face='bold',size=12), 
  axis.text.x=element_text(colour="Black",face="bold",size=10),
  axis.title.x=element_blank(),
  panel.grid.major.x=element_blank(), 
  panel.grid.major.y=element_blank(),
  panel.grid.minor.x=element_blank(), 
  panel.grid.minor.y=element_blank(), 
  panel.background=element_rect(color='black',fill=NA)
) + ylab("log10 (RPM Ecad- / RPM Ecad+") )

snai1 <- assay(ddsMF['Snai2'])
write.csv(snai1, "~/Desktop/Snai1raw.csv")
cdh1_tab$snai1_rpmratio <- scan()
(snai1plot <- ggplot(cdh1_tab, aes(conds,log10(snai1_rpmratio))) +
   ylab("log10 (RPM Ecad- / RPM Ecad+)")+
  geom_boxplot(outlier.color=NA) + ratio_theme +
  ggtitle("Sparc expression\np=0.095") +
  geom_dotplot(binaxis='y', binwidth=0.025, stackdir='center',fill=NA))
  
ggsave(snai2plot, width=3, height=5, file='2014-10-10-snai2-rpm.png')
save(snai1plot, file="2014-10-10-ggplot2_log10targetlevels.RData")
save.image("2014-10-10-batch2.RData")
