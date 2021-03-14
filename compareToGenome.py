# -*- coding: utf-8 -*-
# @Author: Hangyang Zhang
# @Date:   2019-12-09 19:23:48
# @Last Modified by:   Ein
# @Last Modified time: 2019-12-09 19:23:48

fci = 32089.1 # the id of the first chromesome
tata = []; chr_tata = []; current_chr = "1"
with open("v2outWithReverse1000time.gff3") as f:
    for line in f.readlines():
        if line[0] is "#": continue
        ta = line.strip().split()
        if ta[3] is current_chr:
            chr_tata.append(ta)
        else:
            current_chr = ta[3]
            tata.append(chr_tata)
            chr_tata = []; chr_tata.append(ta)
    tata.append(chr_tata)

genes = []
with open("CA_geno.gff") as f:
    for line in f.readlines():
        if "\tgene\t" in line:
            genes.append(line.strip().split())

current_chr = 0; i = 0; re = []; count = len(genes)
for j in range(count):
    overlap_num = 0
    chr = int(float(genes[j][0].split("_")[1]) - fci)
    if chr is not current_chr:
        current_chr = chr; i = 0
    if genes[j][6] is "+":
        while(i < len(tata[chr]) and int(tata[chr][i][4]) < int(genes[j][4])):
            while(i < len(tata[chr]) and tata[chr][i][5] is "+"):
                rl = int(genes[j][3]) - int(tata[chr][i][4]) # relative location
                re.append(["chr"+tata[chr][i][3], tata[chr][i][4], str(rl), tata[chr][i][1],tata[chr][i][2], genes[j][6], genes[j][8].split(";")[2].split("=")[1], tata[chr][i][6]])
                i += 1
            i += 1
    else:
        if j == (count-1): theEnd = 1e08
        else: theEnd = int(genes[j+1][3])
        while(i < len(tata[chr]) and int(tata[chr][i][4]) < theEnd):
            while(i < len(tata[chr]) and int(tata[chr][i][4]) > int(genes[j][3]) and tata[chr][i][5] is "-"):
                rl = int(tata[chr][i][4]) - int(genes[j][4]) # relative location
                re.append(["chr"+tata[chr][i][3], tata[chr][i][4], str(rl+11), tata[chr][i][1],tata[chr][i][2], genes[j][6], genes[j][8].split(";")[2].split("=")[1], tata[chr][i][6]])
                i += 1
                if int(tata[chr][i][4]) > int(genes[j][4]):
                    overlap_num += 1
            if int(tata[chr][i][4]) > int(genes[j][4]):
                overlap_num += 1
            i += 1
        i -= overlap_num

with open("relativeLocation.txt","w") as f:
    f.write("chr\ttataSite\trelativeSite\t-log(score)\tpValue\tchain\tgeneid\ttataseq\n")
    for line in re:
        f.write("\t".join(line)+"\n")
