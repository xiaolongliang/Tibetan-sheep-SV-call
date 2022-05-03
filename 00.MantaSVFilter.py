import sys,os
import re
import gzip

inputs = sys.argv[1]
#name = re.search(r'(.*?)-joint-smoove.genotyped.vcf.gz',inputs)[1]
files = str(inputs) + "/results/variants/diploidSV.vcf.gz"
output = "filter_directory/" + str(inputs) + ".Manta.filter.vcf.gz"
f = gzip.open(output,"wt") 
with gzip.open(files,"rt") as g:
    for line in g:
        line = line.strip()
        if line.startswith("#"):
            f.write(line + "\n")
            continue
        content = line.split()
        if int(content[5]) < 20:
            continue
        f.write(line + "\n")

f.close()
