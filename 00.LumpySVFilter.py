import sys,os
import re
import gzip

inputs = sys.argv[1]
name = re.search(r'(.*?)-joint-smoove.genotyped.vcf.gz',inputs)[1]
output = name + ".genotyped.filter.vcf.gz"
f = gzip.open(output,"wt") 
with gzip.open(inputs,"rt") as g:
    for line in g:
        line = line.strip()
        if line.startswith("#"):
            f.write(line + "\n")
            continue
        content = line.split()
        if content[5] == "0":
            continue
        f.write(line + "\n")

f.close()
