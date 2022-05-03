data <- read.table("Sheep.weir.fst",header=TRUE)
data <- data[!is.na(data$WEIR_AND_COCKERHAM_FST),]
top5 <- quantile(data$WEIR_AND_COCKERHAM_FST,0.95)
fsttop5 <- subset(data,data$WEIR_AND_COCKERHAM_FST > top5)
write.table(fsttop5,file="fsttop0.05.txt",row.names=FALSE,quote=FALSE,sep="\t")
