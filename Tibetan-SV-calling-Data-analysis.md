> Markdown插入图片

![这里写图片描述](https://avatars.githubusercontent.com/u/31619807?s=400&u=0cf97324c487e1c9e3e682427fb81f9d7ab5faa9&v=4)


### 1. 过滤假阳性区域(筛选掉QUAL为LowQD或LowQUAL的SV,只保留PASS的SV)，提取长度和类型
```
python 01.getlength.py vcftools-filter.recode.vcf
```
生成 **Sheep.filter.vcf**
 
### 2. 注释vcf文件（Sheep.filter.vcf）
最后生成 **Sheep_anno.vcf**

### 3. Fst计算
取前5%的Fst结果：
**fsttop0.05.txt**  
前5%对应的基因： 
**GetFstGene.txt**  
对GetFstGene.txt文件进行清理，只保留genes，最后进行富集分析
```
python 00.gene.clean.py GetFstGene.txt | awk '!a[$0]++' | grep -v "LOC" > gene2enrich.txt
```
用于富集分析的纯基因文件：
**gene2enrich.txt**

### 4. SV重复序列注释
INS的fastq序列已经存在在vcf中；  
bedtools提取vcf文件中DEL的fastq
```
bedtools getfasta -fi Texel.fna -bed DEL.bed > DEL.fastq
```
Del.bed文件：chr start end  

用 **RepeatMasker** 软件注释DEL和INS：
```
$ cat DEL.sh
/data/00/software/singularity/singularity_built/bin/singularity exec -B /data -B /home /data/00/software/singularity/TETools RepeatMasker -pa 5 -species human -nolow -no_is -norna -q -dir ./DEL_repeatmasker_repeat DEL.fastq 1>reptmasker.log.o.txt 2>reptmasker.log.e.txt

$ cat INS.sh
/data/00/software/singularity/singularity_built/bin/singularity exec -B /data -B /home /data/00/software/singularity/TETools RepeatMasker -pa 5 -species human -nolow -no_is -norna -q -dir ./INS_repeatmasker_repeat INS.fastq 1>reptmasker.log.o.txt 2>reptmasker.log.e.txt
```
重复序列在 **DEL_repeatmasker_repeat** 和 **INS_repeatmasker_repeat**文件夹中

### 5. 群体分析
###### PCA
```
vcftools --vcf ../Sheep_anno.vcf –plink --out Sheep_anno
plink –file Sheep_anno --pca
Rscript pca.R
```

###### STRUCTURE
```
plink --make-bed --geno 0.05 --maf 0.05 --hwe 0.01 --file Sheep_anno --out ./structure/Sheep_anno

for k in 1 2 3 4 5;do admixture --cv Sheep_anno.bed $k| tee log${k}.out;done

grep -h "CV" log*.out (选择CV最小的，可分的最好！)
CV error (K=1): 0.59358
CV error (K=2): 0.67060 ##
CV error (K=3): 0.76209
CV error (K=4): 0.87673
CV error (K=5): 1.03203

awk 'BEGIN{FS=" ";OFS=","}{print $1,$2,$3,$4,$5}' Sheep_anno.2.Q > Sheep_anno.2.Q.csv

cat Sheep_anno.2.Q.csv| sort -t ","  -k 1n -k 2n -k 3n -k 4n > Sheep_anno.2.Q.sorted.csv

Rscript structure.R  
```
