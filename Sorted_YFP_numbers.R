# R script to make plot of total E and M cells in RNAseq project


samples <- c("PD9210","PD9210","PD8640","PD8640",'PD798','PD798','PD2329','PD2329','PD2342','PD2342','PD1849','PD1849','PD2204','PD2204','PD2523','PD2523','PD345','PD345','PD2412','PD2412')
conds <- c("Ecad+","Ecad-")
conds2 <- rep(conds,10)
conds2
counts <- scan()
df <- data.frame(samples,conds2, counts)
emt <- c("Epithelial",'Mesenchymal')
emt <- rep(emt, 10)
df$emt <- emt
dot_col <- factor(c(rep("red",10),rep('blue',10)))

(rc <- ggplot(df, aes(emt, counts)) + geom_boxplot() + geom_dotplot(binaxis='y', stackdir = 'center', binwidth=15000, fill=NA, col=dot_col))

(rc <- rc + ggtitle("Number of YFP+ Epithelial and Mesenchymal cells\np = 0.1425") +
   ylab("Number of YFP+ cells Ecad+ or Ecad-") +
   scale_y_continuous(labels=comma) + 
   theme(
plot.title = element_text(face='bold'),
axis.text.x=element_text(colour="Black",face='bold',size=15),
axis.text.y=element_text(colour="Black"),
axis.title.x=element_blank(),
axis.title.y=element_text(face='bold'),
panel.grid.major.x=element_blank(),
panel.grid.major.y=element_blank(),
panel.grid.minor.x=element_blank(),
panel.grid.minor.y=element_blank(),
panel.background=element_rect(color='black',fill=NA)
))


ttest <- t.test(Ecad_positive,Ecad_negative)
ttest
setwd("~/Desktop/")
ggsave(rc, file="2014-10-13-number_of_sortedYFP.png")
save.image("2014-10-13-number_of_sortedYFP.RData")

Welch Two Sample t-test

data:  Ecad_positive and Ecad_negative
t = -1.5845, df = 10.593, p-value = 0.1425
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -427949.03   70669.03
sample estimates:
  mean of x mean of y 
141630    320270 