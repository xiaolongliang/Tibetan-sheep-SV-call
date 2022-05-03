#!/bin/bash/python3
import re
import sys,os

"""
usage: python 00.GetFstGene.py ../gene_anno/sv_genes fsttop0.05.txt
"""
def Gene(file1):
    AllGene = {}
    with open(file1,"r") as f:
        for line in f:
            line = line.strip()
            content = line.split()
            if len(content) == 4:
                content.append("NA")
            else:
                pass
            ID = content[2][0:-3]
            #print(content[4])
            TYPE = content[3][1:4]
            AllGene[ID] = [content[0],content[1],TYPE,content[4]]
    return AllGene

def main(file2):
    AllGene = Gene(file1)
    with open(file2,"r") as f:
        for line in f:
            line = line.strip()
            content = line.split()
            tem = content[0] + ":" + content[1]
            if tem in AllGene:
                print("\t".join(AllGene[tem]))

if __name__ == "__main__":
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    main(file2)
