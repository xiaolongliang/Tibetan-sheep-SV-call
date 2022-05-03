#!/usr/bin/python3
#!from ax
"""
Run:
python 01.getlength.py vcftools-filter.recode.vcf
"""
import sys
import re

infile = sys.argv[1]

name = re.search(r'(.*)\.vcf',infile)[1]
outfile = 'Sheep.filter.vcf'
lensize = "Sheep.svlen" + '.txt'

o = open(outfile, 'w')
l = open(lensize, 'w')

with open(infile,'r') as f:

    for line in f:
        line = line.strip()
        if line.startswith('#'):
            o.write(line + '\n')
            continue
        content = line.split('\t')
        if content[6] != "PASS":
            continue
        if re.search(r'SVTYPE=BND', content[7]):
            o.write(line + '\n')
            l.write('BND' + '\t' + '0' + '\n')
            continue
        if re.search(r':SVSIZE=0', content[4]):
        #length = re.search('<(.*):SVSIZE=(\d+)',content[4]).group(2)
            svtype = re.search(';SVTYPE=(.*?);',content[7]).group(1)
            lenend = re.search(';END=(.*?);',content[7]).group(1)
            lenstart = content[1]
            length = int(lenend) - int(lenstart)        
            if length <= 3000 and length > 50:
                o.write(line + '\n')
                l.write(svtype + '\t' + str(length) + '\n')
                continue

        content = line.split('\t')
        length = re.search('<(.*):SVSIZE=(\d+)',content[4]).group(2)
        svtype = re.search('<(.*):SVSIZE=(\d+)',content[4]).group(1)
        length = int(length)
        if length <= 3000 and length > 50:
            o.write(line + '\n')
            l.write(svtype + '\t' + str(length) + '\n')

