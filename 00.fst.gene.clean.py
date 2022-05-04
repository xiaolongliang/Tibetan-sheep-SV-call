#!/bin/python3

"""
To find out the pure gene from 0.05.top.gene 
"""
import re
import sys,os

def main(file):
    gene = {}
    with open(file,"r") as f:
        for line in f:
            line = line.strip()
            content = line.split("\t")
            chr = content[0]
            pos_start = content[1]
            pos_end = content[2]
            body = content[3] + ";"
            gene = re.search(r'(.*?):(.*?);',content[3])
            if gene != None:
                gene = gene.group(2)
                ge = gene.split(",")
                for i in ge:
                    print(i.split("-")[1])
            else:
                continue
            #print("\t".join(gene))

if __name__ == "__main__":
    file = sys.argv[1]
    main(file)
