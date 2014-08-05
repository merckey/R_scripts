options ( stringsAsFactors = FALSE) ;
library ( OmicCircos ) ;
set.seed (1234) ;
## initial
## 10 segments
seg.num <- 10;
## 20 individuals
ind.num <-20; 
## range of the point number is from 20 to 50
seg.po <- c(20:50) ; 
link.num <- 10; 
link.pg.num <-10;
sim.out <- sim.circos ( seg=seg.num , po=seg.po , ind=ind.num , link=link.num , link.pg=link.pg.num);
seg.f <- sim.out$seg.frame ;
seg.v <- sim.out$seg.mapping ;
link.v <-sim.out$seg.link;
link.pg.v <- sim.out$seg.link.pg;
seg.num <- length(unique(seg.f[,1]));
## select segments
seg.name <- paste ("chr" , 1: seg.num , sep="") ; 
db <- segAnglePo ( seg.f , seg=seg.name ) ;
colors <- rainbow(seg.num, alpha=0.5);

