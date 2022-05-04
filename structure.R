png("structure2.png",width=1000,height=600,res=72*2)
dta <- read.csv("Sheep_anno.2.Q.sorted.csv",header=FALSE)
barplot(t(as.matrix(dta)),col=rainbow(2),border=NA)
dev.off()

png("structure3.png",width=1000,height=600,res=72*2)
dta3 <- read.csv("Sheep_anno.3.Q.sorted.csv",header=FALSE)
barplot(t(as.matrix(dta3)),col=rainbow(3),border=NA)

png("structure4.png",width=1000,height=600,res=72*2)
dta4 <- read.csv("Sheep_anno.4.Q.sorted.csv",header=FALSE)
barplot(t(as.matrix(dta4)),col=rainbow(4),border=NA)

png("structure5.png",width=1000,height=600,res=72*2)
dta5 <- read.csv("Sheep_anno.5.Q.sorted.csv",header=FALSE)
barplot(t(as.matrix(dta5)),col=rainbow(5),border=NA)
