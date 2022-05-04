# Data download
本次基于**重测序数据**的SV分析数据是从NCBI上下载的；下载地址为： [SRA数据下载](https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=SRP066883&o=acc_s%3Aa&s=SRR3193935,SRR3193937,SRR3193945)

总共有高低海拔两组羊，高海拔藏羊分为：  
- Prairie Tibetan sheep(PT)  
- Oula sheep(OL)  
- Valley Tibetan sheep(VT)  

低海拔羊分为：  
- Hu sheep(Hu)  
- Small Tail Han sheep(STH);  
    
五种羊每种都下载了10个样本,总共50个样本  

参考基因组为： Texel sheep；下载地址为：[Texel reference genome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/298/735/GCF_000298735.2_Oar_v4.0/)

```shell
## 1.传统方式下载
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR3023075/SRR3023075.1
## 下载文件格式为SRR3023075.1（SRR格式）
## 将SRA转成fastq格式
fastq-dump --split-files SRR3023075.1

## 2.利用ascp高速下载（有时候不能用）
ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR302/005/SRR3023075/SRR3023075_1.fastq.gz .
```

Aspera 下载fastq数据教程： `http://lib.mimazi.net/bioinf/1737.html`

# Align to reference
```shell
bwa mem -t 10 -R "@RG\tID:Hu075\tPL:illumina\tLB:Hu075\tSM:Hu075" ../Texel.fna SRR3023075.1_1.fastq SRR3023075.1_2.fastq | samtools view -Sb - > ./BAM_directory/Hu075.bam
samtools sort Hu25.bam -o Hu25.sorted.bam && samtools index Hu25.sorted.bam
```
# CALL SV
## 1. Delly

```shell
## download Delly by bioconda
conda install -c bioconda delly
## call sv for every individual
delly call -g ../Texel.fna -o ./Delly_directory/Hu075.bcf ./BAM_directory/Hu075.sorted.bam
## merge sv for all individual
delly merge -o all.sites.bcf `cat 00.bcflist.oneline`
## call sv in every site for every samples
delly call -g ../../Texel.fna -v all.sites.bcf -o Hu075.geno.bcf ../BAM_directory/Hu075.sorted.bam
## convert geno.bcf to geno.vcf
bcftools view Hu075.geno.bcf -O v > Hu075.vcf
## filter sv that IMPRECISE and LowQual
python 00.DellySVFilter
## the final sv file for delly call is Hu075.filter.vcf 
```

>cat 00.bcflist.oneline  
>Hu075.bcf Hu25.bcf Hu52.bcf Hu53.bcf Hu75.bcf Hu76.bcf Hu77.bcf Hu84.bcf Hu94.bcf... 
## 2. Lumpy
```shell
## download Lumpy(depend on python2.7)
conda install -c bioconda lumpy-sv
## activite the Lumpy environment
conda activate /data/00/user/user211/miniconda3/envs/PYTHON2
## call sv for every individual
smoove call --outdir ./Lumpy_directory --name Hu075 --fasta ../Texel.fna -p 10 --genotype ./BAM_directory/Hu075.sorted.bam
## combined all vcf
smoove merge --name merged -fasta ../../Texel.fna *.genotyped.vcf.gz
## call sv in every sites for every individual
smoove genotype -d -x -p 10 --name Hu075-joint --outdir results-smoove/ --fasta ../../Texel.fna --vcf merged.sites.vcf.gz ../BAM_directory/Hu075.sorted.bam
## filter out QUAL=0
python 00.LumpySVFilter
```
## 3. Manta
```shell
python /data/01/user196/work2021/sheep/BIG_89sheep/manta-1.6.0.centos6_x86_64/bin/configManta.py --bam Hu075.sorted.bam --referenceFasta ../Texel.fna --runDir ./Manta_directory/Hu075
python ./Manta_directory/Hu075/runWorkflow.py
## filter out the sv
sh 00.MantaSVFilter.run.sh
```
# Combined the sv results worked by three software
1. svimmer combined  
```shell
## bgzip压缩，tabix建立索引
for i in *.vcf;do bgzip $i;done
for i in *.vcf.gz;do tabix -p vcf $i;done
svimmer --threads 20 all.vcf.list LG01 LG02 LG03 LG04 ......> input.vcf
bgzip input.vcf
tabix -p vcf input.vcf.gz
## all.vcf.list中包括三个软件所有的过滤后的vcf结果文件
## 染色体名称必须在参考基因组的vcf文件中提取（awk）
```
2. graphtype SV genotype again
```shell
## 以NC_019458.2这条染色体为例
graphtyper genotype_sv ../../Texel.fna input.vcf.gz --sams=bam.list --region=NC_019458.2 1 275406953
## bam.list中包括所有需要合并的个体的bam文件
## 脚本中使用的chrpos.bed需要自己生成, chrpos.bed文件的每一行包含每条染色体的名称以及起始位置
awk '{print $1, 1, $2}' cishu.LG.fasta.fai > chrpos.bed
## 最后自动生成文件夹sv_results，所有的文件都在此文件夹下
```
```python
## 00.graphtyper.py脚本可批量生成每条染色体的运行命令
import os
import re
with open("chrpos.bed","r") as f:
    for line in f:
        line = line.strip()
        print("graphtyper genotype_sv ../../Texel.fna input.vcf.gz --sams=bam.list --region=" + line)
```
3. combined vcf for every chromsome
```shell
cd sv_results
grep '^#' chr1.vcf > merge.vcf
grep -v '^#' chr1.vcf  chr2.vcf chr3.vcf chr4.vcf >> merged.vcf
```
对于此项目，因为基因组中contig很多，所以会生成很多的contig对应的vcf文件， 写shell脚本提取所有的染色体和contig对应的vcf
```shell
#!/bin/bash
for i in `ls`;do
for j in ${i}/*vcf.gz;do
if [ -f "$j" ];then  # -f:判断文件是否存在
echo $j >> vcf.list
else
continue
fi
done
done
```
4. 使用bcftools去掉merged.vcf中sv的BREAKPOINT
```shell
bcftools view --include "SVMODEL='AGGREGATED'" -Oz -o get_agg.vcf.gz merged.vcf
```
5. 使用vcftools过滤掉GT的max missing0.95的SV
```shell
vcftools --vcf get_agg.vcf.gz --max-missing 0.95 --recode --recode-INFO-all --out vcftools-filter
```
6. 过滤假阳性区域，提取长度和类型
```python
python 00.filterSV&getTypeLen vcftools-filter.recode.vcf
## 最后生成文件： Sheep.filter.vcf
```
# sv annotation
sv 注释用到的软件为 [vcfanno](https://github.com/brentp/vcfanno)
```shell
# 1. 生成文件
perl -alne 'if ($F[2] eq "gene") {$F[8] =~ /ID=(.*?);/; $name = $1; $start = $F[4] - 1; $end = $start + 20001; print "$F[0]\t$start\t$end\t$name";}' GCF_000298735.2_Oar_v4.0_genomic.gff | sort -k1,1 -k2,2n -k3,3n | bgzip -c > Texel.genes.downstream.bed.gz
perl -alne 'if ($F[2] eq "gene") {$F[8] =~ /ID=(.*?);/; $name = $1; $start = $F[3] - 20001; $end = $F[3]; print "$F[0]\t$start\t$end\t$name";}' GCF_000298735.2_Oar_v4.0_genomic.gff | sort -k1,1 -k2,2n -k3,3n | bgzip -c > Texel.genes.upstream.bed.gz
perl -alne 'if ($F[2] eq "exon") {$F[8] =~ /ID=(.*?);/; $name = $1; $start = $F[3] - 1; $end = $F[4] - 1; print "$F[0]\t$start\t$end\t$name";}' GCF_000298735.2_Oar_v4.0_genomic.gff | sort -k1,1 -k2,2n -k3,3n | bgzip -c > Texel.exons.srt.bed.gz
perl -alne 'if ($F[2] eq "gene") {$F[8] =~ /ID=(.*?);/; $name = $1; $start = $F[3] - 1; $end = $F[4] - 1; print "$F[0]\t$start\t$end\t$name";}' GCF_000298735.2_Oar_v4.0_genomic.gff | sort -k1,1 -k2,2n -k3,3n | bgzip -c > Texel.genes.srt.bed.gz
perl -alne 'if ($F[2] eq "CDS") {$F[8] =~ /ID=(.*?);/; $name = $1; $start = $F[3]; $end = $F[4]; print "$F[0]\t$start\t$end\t$name";}' GCA_017524585.1_CAU_O.aries_1.0_genomic.gff | sort -k1,1 -k2,2n -k3,3n | bgzip -c > Tibetan.cds.srt.bed.gz

python cds.py  | bgzip -c  > Tibetan.cds_end.srt.bed.gz
python cds.py  | bgzip -c  > Tibetan.cds_start.srt.bed.gz

# 2. 建立索引
tabix Tibetan.genes.upstream.bed.gz
tabix Tibetan.genes.downstream.bed.gz
tabix Tibetan.cds_end.srt.bed.gz
tabix Tibetan.cds_start.srt.bed.gz
tabix Tibetan.exons.srt.bed.gz
# 3. 注释
# sv number
./vcfanno_linux64 config.toml ../Sheep.filter.vcf > Sheep_anno.vcf
```

# Fst
```shell
## calculate the fst between low altitude sheep and Tibetan sheep
vcftools --vcf ../gene_anno/Sheep_anno.vcf --weir-fst-pop low_sheep --weir-fst-pop high_Tibetan --out Sheep
## select the top 5% 
Rscript fstTop0.05.R
## the top 5% fst was saved by fsttop0.05.txt
## get the genes based on fst top 5%
python 00.GetFstGene.py ../gene_anno/sv_genes fsttop0.05.txt  > GetFstGene.txt
```  
[群体分析详细过程](https://github.com/xiaolongliang/Tibetan-sheep-SV-call/blob/master/Tibetan-SV-calling-Data-analysis.md) 

###### 补充：
```shell
## 从注释后的vcf文件Sheep_anno.vcf中中提取每个SV对应的基因
python SV_gene.py > sv_genes
## 从sv_genes中提取每个SV在基因上的位置，包括downstream, upstream, genes, exons等
python 00.Getgenetype.py > Getgenetype
```
