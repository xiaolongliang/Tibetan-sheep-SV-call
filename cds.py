#!/bin/python3

import gzip
import re,sys

def main(file):
    cds = {}
    with gzip.open(file,"rt",encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            content = line.split()
            if content[3] not in cds:
                cds[content[3]] = [[content[1],content[2],content[0]]]
            else:
                cds[content[3]].append([content[1],content[2],content[0]])
    return cds

def run():
    cds = main(file)
    cds_out = {}
    for i in cds:
        if len(cds[i]) == 1:
            cds_start = cds[i][0][0]
            cds_end = cds[i][0][1]
            cds_id = cds[i][0][2]
        else:
            cds_start = cds[i][0][0]
            cds_end = cds[i].pop()[1]
            cds_id = cds[i][0][2]

        #generate the Tibetan.cds.srt.bed.gz   "python cds.py  | bgzip -c  > Tibetan.cds_start.srt.bed.gz" 
        print("\t".join([cds_id,str(int(cds_start)-1),str(int(cds_start)),i]))
        #generate the Tibetan.cds_end.srt.bed.gz   "python cds.py | bgzip -c  > Tibetan.cds_end.srt.bed.gz"
        #print("\t".join([cds_id,str(int(cds_end)-1),str(int(cds_end)),i]))  

if __name__ == "__main__":
    file = "Texel.cds.srt.bed.gz"
    run()
