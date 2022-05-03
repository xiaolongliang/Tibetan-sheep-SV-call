import sys,re
infile = sys.argv[1]
name = re.search(r'(.*)\.vcf',infile)[1]
outfile = name + ".filter.vcf"
o = open(outfile,"w")
with open(infile,"r") as f:
    for line in f:
        line = line.strip()
        if line.startswith("#"):
            o.write(line + "\n")
            continue
        content = line.split()
        if content[6] != "PASS" or re.match(r"IMPRECISE",content[7]):
            continue
        o.write(line + "\n")
o.close()
