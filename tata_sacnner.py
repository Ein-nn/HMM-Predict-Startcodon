# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 08:38:33 2019
@author: Hanyang Zhang
@ID: 1730416017
"""
import random
import math
class TataScaner:
    def __init__(self, HMMmodel = {}, seqs = [], cutoff = 0.05, times = 100):
        self.hmm = HMMmodel
        self.seqs = seqs
        self.cutoff = cutoff
        self.score = [] 
        self.times = times 
    
    def shuffleStr(self, seq=""):
        seqList = list(seq)
        random.shuffle(seqList)
        return "".join(seqList)
    
    def reverse(self, seq=""):
        dic = {'A':'T','T':'A','C':'G','G':'C'}
        re = ''
        for i in range(len(seq)):
            re += dic[seq[len(seq)-i-1]]
        return re
    
    def loadHmmFromFile(self, filename = ""):
        with open(filename,"r") as f:
            for line in f.readlines():
                if line[0] is '%':
                    data = line.split()
                    for i in range(1,13):
                        self.hmm[str(i)+data[0].replace('%','')] = round(float(data[i])/100, 3)
    
    def loadSeqsFromFile(self, filename = ""):
        seq = ""; self.chrome = []
        with open(filename, "r") as f:
            for line in f.readlines():
                if line[0] is ">":
                    self.seqs.append(seq)
                    chr = line.replace("Candida albicans SC5314 chromosome","").replace("sequence","").strip()
                    self.chrome.append(chr[1:].split()); seq = ""
                else:
                    seq += line.strip().upper()
            self.seqs.append(seq); self.seqs.remove("")
        
    def countScore(self, seq = ""):
        s = 1
        for i in range(len(seq)):
            s *= self.hmm[str(i+1)+seq[i]]
        return s
    
    def stepScore(self, seq = "", chr = 0, site = 1):
        composition = seq.count('A')+seq.count('T')+seq.count('C')+seq.count('G')
        if composition != len(seq):
            return None
        s = self.countScore(seq) 
        pValue = self.countP(seq, s)
        if pValue is not None:
            if pValue < self.cutoff:
                self.score.append([str(s),str(-math.log10(s)),str(pValue),str(self.chrome[chr][1]),str(site),"+",seq])
        
        seq = self.reverse(seq) 
        s = self.countScore(seq)
        pValue = self.countP(seq, s)
        if pValue is not None:
            if pValue < self.cutoff:
                self.score.append([str(s),str(-math.log10(s)),str(pValue),str(self.chrome[chr][1]),str(site),"-",seq])
    
    def countP(self, seq = "", s = 0): 
        if s == 0:
            return None
        if self.repeatCount(seq, s, 10) > 0.5: 
            return None
        return self.repeatCount(seq, s, countTime = self.times)
    def repeatCount(self, seq = "", s = 0, countTime = 100): 
        randomSeqs = {}; n = 0; maxTime = countTime * 2; time = 0
        while(len(randomSeqs.keys()) < countTime and time < maxTime):
            randomseq = self.shuffleStr(seq)
            time += 1
            if randomseq in randomSeqs:
                continue
            randomScore = self.countScore(randomseq)
            randomSeqs[randomseq] = randomScore
            if self.countScore(randomseq) > s:
                n += 1
        return (n/len(randomSeqs.keys()))
    
    def scan(self, wordSize = 12, step = 1): 
        for i in range(len(self.seqs)):
            for j in range(len(self.seqs[i]) - wordSize + step):
                self.stepScore(self.seqs[i][j:j+wordSize], chr = i, site = j+1)
def run():
    ts = TataScaner(cutoff = 0.005, times = 1000)
    ts.loadHmmFromFile('./TATA-boxHMM.txt')
    ts.loadSeqsFromFile('./test.fna')
    ts.scan(wordSize=12, step=1)
    with open("test.gff","w") as f:
        f.write("#Score\t-log(score)\tp-value\tchr\tstartSite\tchain\tseq\n")
        for i in ts.score:
            f.write("\t".join(i)+"\n")
if __name__ == '__main__':
    run()
