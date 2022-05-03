#!/bin/python3

import sys,os
import re

file = "sv_genes"
with open(file,"r") as f:
    for line in f:
        line = line.strip().split()
        if len(line) != 5:
            continue
        body = line[4] + ";"
        Getype = re.findall(r"(.*?):(.*?);",body)
        Getype = [i[0] for i in Getype]
        for j in Getype:
            print(line[0] + "\t" + line[1] + "\t" + j)
