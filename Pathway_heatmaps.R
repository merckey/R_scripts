
library(ggplot2)
library(RColorBrewer)
library(reshape)


df <- read.delim(file = "2014-10-13-DESeq2-ddsMF_8_Axon.csv", header=T)
head(df)

# generate Ecad-/Ecad+ ratios 
mat <- (df$PD2204E_minus/df$PD2204E_plus)
mat2 <- (df$PD2329E_minus/df$PD2329E_plus))
mat2 <- (df$PD2329E_minus/df$PD2329E_plus)
mat3 <- (df$PD2342E_minus/df$PD2342E_plus)
mat4 <- (df$PD2412E_minus/df$PD2412E_plus)
mat5 <- (df$PD2523E_minus/df$PD2523E_plus)
mat6 <- (df$PD345E_minus/df$PD345E_plus)
mat7 <- (df$PD798E_minus/df$PD798E_plus)
mat8 <- (df$PD9210E_minus/df$PD9210E_plus)

df <- data.frame(df$Gene,mat, mat2,mat3,mat4,mat5,mat6,mat7,mat8)
colnames(df) <- c('Gene',"PD2204","PD2329",'PD2342','PD2412',
                   'PD2523','PD345','PD798','PD9210')

head(df)

dfMelt <- melt(df)
head(dfMelt)

# heatmaps
(q <- ggplot(dfMelt, aes(variable, Gene, fill=log(value,2))) + geom_tile(color='black') +
  scale_fill_gradient2(low='blue',mid="white",high='red',midpoint=0))
(q + theme_bw())

ggsave(q, "2014-10-14-Axon-pathway-heatmap.png")
save.image("2014-10-14-Axon-pathway-heatmap.RData")
