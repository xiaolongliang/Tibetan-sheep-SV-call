#!/bin/python3

"""
The script is used to get the genes(included the exon,upstream,downstream...) of 
each sv from sv annotated file named Sheep_anno.vcf(annotated by config.toml; get by anxuan);
"""

import re
import sys,os

def main(file):
    with open(file,"r") as f:
        for line in f:
            line = line.strip()
            content = line.split()
            if line.startswith("#"):
                continue
            body = content[7] + ";"
            CHR = content[0]
            POS = content[1]
            ID = content[2]
            ALT = content[4]
            GENE = re.findall(r"Texel_(.*?)=(.*?);",body)
            gene = [":".join(i) for i in GENE]
            print("\t".join((CHR,POS,ID,ALT,";".join(gene))))

if __name__ == "__main__":
    file = "Sheep_anno.vcf"
    main(file)
