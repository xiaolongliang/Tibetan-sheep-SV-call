rm(list = ls())
getwd()
setwd("E:/#课题/藏绵羊/结果/20220318_50samples")

# 1. plot the sv number
library(tidyverse)
data <- read.table("svlen.bar")
colnames(data) <- c("Svtype","Svlen")
data %>% add_row(Svtype = "Total",Svlen=sum(data$Svlen)) -> data

data$Svtype <- factor(data$Svtype,levels = c("DEL","INS","BND","INV","DUP","Total"))

# theme_set(theme_classic())
ggplot(data = data,mapping = aes(x=reorder(Svtype,-Svlen),y=Svlen,fill=Svtype)) + 
  geom_bar(stat = "identity",position = "dodge") + 
  coord_flip() +
  labs(x = "",y="SV number") + 
  theme_bw() 
ggsave(filename = "SVtypenumber.tiff",height = 6,width = 9,dpi = 300)  


# 2. sv filter based on frequence
freq <- read.table("sv_freqs.txt",header = TRUE)

freq[freq=="0/0"] = as.numeric(0)
freq[freq=="1/0"] = as.numeric(1)
freq[freq=="0/1"] = as.numeric(1)
freq[freq=="1/1"] = as.numeric(2)
freq[freq=="./."] = NA

Tibetan <- freq %>% select(contains(c("PT","Oula","VT"))) 
Tibetan <- as.data.frame(lapply(Tibetan, as.numeric)) #dataframe中的字符串转变成数字
lowSheep <- freq %>% select(contains(c("Hu"))) #"STH"
lowSheep <- as.data.frame(lapply(lowSheep, as.numeric)) 

# subset(freq,TibetanF >= 0.7 & LowsheepF <= 0.3,select = c(CHROM,POS,TibetanF,LowsheepF))
# prop.table(table(Tibetan$PT35)) #table();   prop.table(table(x)):计算频率 

# calculate the frequence
g <- function(x){
  out <- c()
  for(i in 1:nrow(x)){
    f = prop.table(table(t(x[i,])))
    tem = f[names(f) == "2"]
    if(length(tem) != 0){
      out <- c(out,tem)
    }else{
      out <- c(out,0)
    }
  }
  return(out)
}

Tibetanfreq_2 <- g(Tibetan)
Tibetanfreq_2 <- as.data.frame(Tibetanfreq_2)
Lowsheepfreq_2 <- g(lowSheep)
Lowsheepfreq_2 <- as.data.frame(Lowsheepfreq_2)

FreqFilterOut <- cbind(freq,Tibetanfreq_2,Lowsheepfreq_2)
# result1 <- subset(FreqFilterOut,Tibetanfreq_2 >= 0.3 & Lowsheepfreq_2 <= 0.7,select = c(CHROM,POS,Tibetanfreq_2,Lowsheepfreq_2,genes))
# result2 <- subset(FreqFilterOut,Tibetanfreq_2 <= 0.7 & Lowsheepfreq_2 >= 0.3,select = c(CHROM,POS,Tibetanfreq_2,Lowsheepfreq_2,genes))
# result <- rbind(result1,result2)
write.table(result,file = "SvFilter2Freq",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(FreqFilterOut,file = "SvFilter2FreqAll.txt",quote = FALSE,sep = "\t",row.names = FALSE)
# aa <- FreqFilterOut %>% select(CHROM,POS,Tibetanfreq_2,Lowsheepfreq_2)

result <- subset(FreqFilterOut,abs(Tibetanfreq_2 - Lowsheepfreq_2) >= 0.7,select = c(CHROM,POS,Tibetanfreq_2,Lowsheepfreq_2,genes))

# 3. plot the sv repeat
data <- read.csv("sv_repeat.csv",header = FALSE,encoding = "UTF-8")
names(data) <- c("RepeatType","number","svType")

theme_set(theme_classic())
ggplot(data,aes(x=svType,y=number,fill=RepeatType)) + geom_bar(position = "stack",stat = "identity") +
  labs(y = "Repeat Number",x = "SV Type",fill = "Repeat Type")

ggsave("SV_Repeat.tiff",height = 6,width = 8,dpi = 200)



## PANTHER
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("PANTHER.db")

# 4. sv length plot 
len <- read.table("Sheep.svlen.txt",header = FALSE)

# len %>% filter(V1=="DEL") %>% group_by(V2) %>% summarise(count=n()) -> DEL
# len %>% filter(V1=="INS") %>% group_by(V2) %>% summarise(count=n()) -> INS
# len %>% filter(V1=="BND") %>% group_by(V2) %>% summarise(count=n()) -> BND
# len %>% filter(V1=="INV") %>% group_by(V2) %>% summarise(count=n()) -> INV
# len %>% filter(V1=="DUP") %>% group_by(V2) %>% summarise(count=n()) -> DUP

ggplot(len,aes(x=V2,fill=V1)) + geom_density(alpha=.3) + 
  scale_x_log10() +
  scale_y_continuous(expand = c(0,0),limits=c(0,4)) + 
  labs(x = "Variant Length (bp)",fill="SV Type")
ggsave("sv_length.tiff",height = 6,width = 13,dpi = 150)


# 5. sv gene type
rm(list = ls())

type <- read.table("Getgenetype",header = FALSE)
names(type) <- c("ID","pos","type")
type %>% group_by(type) %>% summarise(n = n()) -> plot
plot %>% add_row(type = "CDS",n = sum(122,128)) -> plot
plot <- plot[-c(1,2),]

ggplot(plot,aes(x = type,y=n,fill=type)) + geom_bar(stat = "identity") + 
  labs(x = "", y = "SV Count")
ggsave("Getgenetype.tiff",height = 10,width = 10,dpi = 150)



# 6. sv PCA
library(ggplot2)
dta <- read.table("plink.eigenvec",header=FALSE,sep=" ")
dta <- dta[1:12]
names(dta) <- c("ID1","ID2",paste0("PC",1:10))
dta$FILL <- c(rep(1,10),rep(2,10),rep(2,10),rep(1,10),rep(2,10))
dta$FILL <- factor(dta$FILL)
#write.table(dta,file="pcafile.txt",quote=FALSE,row.names=TRUE,sep="\t")
b <- ggplot(dta,aes(x = PC1,y = PC2,color=FILL,shape=FILL))
b + geom_point() +
  labs(x = "PC1 (2.37%)", y = "PC2 (2.35%)") +
  theme_light() +
  scale_color_discrete(name="Species",breaks=c(1,2),labels=c("Low Sheeps","Tibetan Sheeps")) +
  scale_shape_discrete(name="Species",breaks=c(1,2),labels=c("Low Sheeps","Tibetan Sheeps"))  
  # scale_y_continuous(expand = c(0,0)) + 
  # scale_x_continuous(expand = c(0,0)) 
  
ggsave("pcaBsv.tiff",heigh=5,width=7,dpi=150)


# 7. FST
library(CMplot)
fst <- read.table("Sheep.weir.fst",header = TRUE)
end <- read.table("sv_end.txt")
names(end) <- "END"
fst <- cbind(fst,end)


fst %>% filter(!is.na(WEIR_AND_COCKERHAM_FST)) -> fst
fst$WEIR_AND_COCKERHAM_FST[fst$WEIR_AND_COCKERHAM_FST < 0] = 0
fst %>% select("CHROM", "POS","END","WEIR_AND_COCKERHAM_FST") -> fst

fst0.1 <-  quantile(fst$WEIR_AND_COCKERHAM_FST,0.99)
unname(fst0.1)


# grep ">" GCF_000298735.2_Oar_v4.0_genomic.fna | head -n 27 | sed -r 's/^>(.*)Ovis.*Texel (.*), Oar_v4.0.*/\1 \2/' > NC2chr
chr <- read.table("NC2chr")
chr$Chr <- paste0(chr$V2,chr$V3)

for (i in 1:nrow(fst)) {
  if (fst$CHROM[i] %in% chr$V1){
    index = match(fst$CHROM[i],chr$V1)
    fst$CHROM[i] = chr$V3[index]
  }else{
    next
  }
}

fst <- data.frame(ID=paste(fst$CHROM,fst$POS,sep = ":"),
                  chr = fst$CHROM,
                  pos = fst$POS,
                  val = fst$WEIR_AND_COCKERHAM_FST)

fst %>% slice(2:47948) -> fst

CMplot(fst,
       type = "p",  
       plot.type = "m",
       LOG10 = F,
       col = "grey",
       cex = 0.2,
       band = 0.5,
       ylab = "Fst",
       ylab.pos = 2.5,
       #ylim = c(0,10),
       cex.axis = 0.8,
       threshold = fst0.1, ## 显著性阈值
       threshold.col = 'red',  ## 阈值线颜色
       threshold.lty = 2, ## 阈值线线型
       threshold.lwd = 2, ## 阈值线粗细
       amplify = F,  ## 放大
       #ylab=expr(paste(pi, "_pop1/", pi, "_pop2")),
       highlight=snp,
       # highlight.col=c("red","blue","green"),
       # highlight.cex=1,
       # highlight.pch=c(15:17), 
       highlight.text=genes,      
       # highlight.text.col=c("red","blue","green"),
       file.output = F)


snp <- c("22:20536060","15:53315375","5:80272460","20:17310296")
genes <- c("HIF1AN","UVRAG","XRCC4","VEGF165")




## 8. Tibetan sheep VS hu sheep 
## blood parameters
setwd("E:/#课题/藏绵羊/结果/藏羊和湖羊的血液参数分析")
getwd() 
dta <- read.table("bloodparameters_Tibetan_vs_husheep.txt",header = TRUE,na.strings = "NA",sep = "\t")
names(dta)[1] <- "Para"

## 单因素方差分析 one-way analysis of variance

RBC <- dta[1,]  |> gather(key = key,value = value,-Para) |> 
  mutate(key = factor(key,levels = c( "WM.maqu","WM.luqu","GJ","SK","Hu"))) -> RBC

aov(value~key,data = RBC) -> aov.manu
summary(aov.manu)

library(multcomp, quietly = TRUE)
glht(aov.manu,linfct = mcp(key = "Tukey")) |> summary()
