#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
import shutil
import gzip
import numpy as np
import time

t1 = time.time()
#import math
#import multiprocessing as mp

#def main():
def count_lines(filename, chunk_size=1<<13):
    with open(filename) as file:
        lenmax=len(max(filename, key = len))
        return sum(chunk.count('\n')
             for chunk in iter(lambda: file.read(chunk_size), ''))
   
    
filename = "SRR11700309.fastq"
file = os.path.join(os.getcwd(), filename)
Nlines=(count_lines(filename),)[0]

with open(file, "r") as infile1:
    lenmax=(len(max(infile1, key = len)))-1    

with open(file, "r") as infile:
#    

#    lenmax=len(max(filename, key = len))
    
    line = infile.readline().strip()
    indseq = []
    linecount = 0
    readcount = 0
    lenmin = 300
    sumlen = 0

    lenmax1 = 0
    lenmean_sum = 0
    array_phreds = np.zeros((Nlines // 4, lenmax))
    array_seq = np.zeros((Nlines // 4, lenmax))
    while linecount < Nlines:

        if linecount % 4 == 0: #Name
            pass
        elif linecount % 4 == 1: #Seq
            seq = [x for x in line]
            seq = []
            for x in line:
                if x == 'A':
                    seq.append(1)
                elif x == 'C':
                    seq.append(2)
                elif x == 'G':
                    seq.append(3)
                elif x == 'T':
                    seq.append(4)
                else:
                    seq.append(5)
            arr_seq = np.array(seq).astype(np.float16)
            Nnans = lenmax - len(arr_seq) 
            
            arr_seq = np.pad(arr_seq, [(0, Nnans)], mode='constant', constant_values=np.nan)

            array_seq[readcount-1] = arr_seq

        elif linecount % 4 == 2: #+
            pass
        elif linecount % 4 == 3: #Phred
            phredstr = [ord(x) - 33 for x in line]
            arr = np.array(phredstr).astype(np.float16)
            Nnans = lenmax - len(arr) 
            arr = np.pad(arr, [(0, Nnans)], mode='constant', constant_values=np.nan)

            array_phreds[readcount] = arr
            if len(phredstr) < lenmin:
                lenmin = len(phredstr)
            sumlen += len(phredstr)
            readcount += 1
        line = infile.readline().strip()
        linecount += 1


print("""

       ***********************************
       *                                 *                                
       *        BASIC STATISTICS         *                      
       *                                 *                              
       ***********************************
""",)
   

print("Filename = {0:s}".format(filename))

max_phred=np.amax(array_phreds[1])
if max_phred>35:
    Sequencing_Platform = "Illumina"
    print("Sequencing Platform = Illumina")
else:
    Sequencing_Platform = "Ion Torrent"
    print("Sequencing Platform = Ion Torrent")

if Sequencing_Platform == "Illumina" or "Ion Torrent":
    encoding = 'Illumina 1.9 / Sanger / Ion Torrent'
    phred = 'Phred+33'    

print("Encoding = {0:s} \nPhred = {1:s}".format(encoding, phred))
print("Total Sequences = {0:d}".format(readcount))

phredcounter = 0
poorseq = 0
while phredcounter < readcount:
    mean_phred=np.mean(array_phreds[phredcounter])
    if mean_phred < 10:
        poorseq += 1
    phredcounter += 1
print('Sequences flagged as poor quality = {0:d}'.format(poorseq))

lenmean = sumlen // readcount
print('Sequence length = {0:d} - {1:d} (mean: {2:d})'.format(lenmin,lenmax,lenmean))

s2=np.sum(array_seq==2)
s3=np.sum(array_seq==3)
sv=array_seq.size
s23=s2+s3
proc= (100*s23)//sv
print('GC% = {0:d}'.format(proc))
print("\n")


t2 = time.time()
print('Working time = {0:f}'.format(t2 - t1))
