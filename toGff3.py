# -*- coding: utf-8 -*-
# @Author: Hangyang Zhang
# @Date:   2019-12-19 08:23:31
# @Last Modified by:   Ein
# @Last Modified time: 2019-12-19 08:23:33

n = 1; cn = "chr1"; chr_id = 32089.1;
with open("TATAfiltered.txt","r") as f:
    for line in f.readlines():
        with open("TATAfiltered.gff3", "a") as w:
            dat = line.split()
            if(dat[0] != cn):
                cn = dat[0] ; chr_id += 1
            chr = "NC_0" + str(chr_id)
            id = "id=TATA-box_"+str(n)+";name=TATA_"+dat[6]
            outdat = [chr,".","TATA-box",dat[1],str(int(dat[1])+11),dat[4],dat[5],".",id]
            w.write("\t".join(outdat)+"\n"); n += 1
