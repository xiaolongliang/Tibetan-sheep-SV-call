import os
import re
with open("chrpos.bed","r") as f:
    for line in f:
        line = line.strip()
        print("graphtyper genotype_sv ../../Texel.fna input.vcf.gz --sams=bam.list --region=" + line)

