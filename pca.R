library(ggplot2)
dta <- read.table("plink.eigenvec",header=FALSE,sep=" ")
dta <- dta[1:12]
names(dta) <- c("ID1","ID2",paste0("PC",1:10))
dta$FILL <- c(rep(1,10),rep(2,10),rep(2,10),rep(1,10),rep(2,10))
dta$FILL <- factor(dta$FILL)
#write.table(dta,file="pcafile.txt",quote=FALSE,row.names=TRUE,sep="\t")
b <- ggplot(dta,aes(x = PC1,y = PC2,color = FILL,shape=FILL))
b + geom_point(size=1) +
	theme_light() +
	#guides(color="none") + 
	scale_color_discrete(name="Species",breaks=c(1,2),labels=c("Low Sheeps","Tibetan Sheeps")) +
	scale_shape_discrete(name="Species",breaks=c(1,2),labels=c("Low Sheeps","Tibetan Sheeps")) +
	labs(x = "PC1 (2.37%)", y = "PC2 (2.35%)")
ggsave("pca.tiff",heigh=7,width=9,dpi=150)
