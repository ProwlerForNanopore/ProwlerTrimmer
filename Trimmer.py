#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 00:50:03 2020

@author: simon
"""
import Fastq as fq
#import pylab as pl
#import matplotlib.pyplot as plt
import time as tm
#import numpy as np
#import collections as cl
import sys

#%%
def readFile(file):
  
  #file = "SRR9987849_1.fastq"
  print("number of seqs: ", fq.getFileLen(file)/4)
  #seqsToAnalyze = int(fq.getFileLen(file)/4)
  
  #reads the fastQ file and creates a list of fastq sequences
  readStart = tm.time()
  
  seqs = fq.readFastqFile(file,seqsToAnalyze)
  print("read time: ", tm.time()-readStart)
  
  return seqs

#%%
def rollingTrim(seqs, mode, Qcutoff, windowSize, minLen):
  # uses rolling window trimmer for each sequence
  # then sets up the rolling quality plot for each sequence
  print("mode = ", mode)
  startTime = tm.time()
  fq.rollingTrimAllSeqs(seqs,Qcutoff,windowSize,mode,minLen)
  endTime = tm.time()
  print("end trim. Time elapsed: ", endTime-startTime)
  
  return seqs

#%%

#contigsSingle = seqs[0].generateContigs(0,0)
#contigsMulti = fq.generateContigsList(seqs[0:seqsToAnalyze],0,0)
#contigLenList = []
#for i in range(len(contigsMulti)):
#  contigLenList.append(len(contigsMulti[i]))
#lenListLib = cl.Counter(contigLenList)
#print(lenListLib)

#%%

def writeFile(file, seqs, mode, Qcutoff, windowSize, minLen, seqsToAnalyze, outMode="fq"):
  # add code suffixes to output file name for mode(D__/S__), window size (W__),
  # minimum Length of read (L__), and sequences selected from raw file (R[__-__)
  startTime = tm.time()
  outputFileName = file[0:-6] + "Trim"
  outputFileName += mode + str(Qcutoff)
  outputFileName += "W" 
  outputFileName += str(windowSize) 
  outputFileName += "L" + str(minLen)
  outputFileName += "R" + str(seqsToAnalyze)
  
  if outMode == "fq":
    fq.writeFastq(seqs,outputFileName)
  elif outMode == "fa":
    fq.writeFasta(seqs,outputFileName)
  else:
    print("error. output file type not supported")
  print("write took: ", tm.time() - startTime)  
  #fq.writeFasta(seqs)

#%%
#for Qcutoff in range(18,19,2):
print("Qcutoff: ",Qcutoff)
if __name__ == "__main__":
#  print(f"Arguments count: {len(sys.argv)}")
#  for i, arg in enumerate(sys.argv):
#      print(f"Argument {i:>6}: {arg}")

  print("program start:")      

  
  #defines the input variables from commandlince arguments
#  file = sys.argv[1]
#  windowSize = int(sys.argv[2])
#  minLen = int(sys.argv[3])
#  mode = sys.argv[4]
#  Qcutoff = int(sys.argv[5])
#  seqsToAnalyze = int(sys.argv[6])
#  outMode = sys.argv[7]

#    file = "G:/Honours Thesis Data/SequenceData/testSeq"
  file = "G:/Honours Thesis Data/SequenceData/testSeq674"
  fileIn = file +".fastq"
  fileOut = file + "Out" + ".fastq"
#    file="G:/Honours Thesis Data/SequenceData/SAMN08158127MesoplasmaSyrphidaeStrainYJS/SRR6364637_1-MesoplasmaSyrphidaeStrainYJS-MinIONRaw.fastq"
  windowSize = 10
  minLen = 10
  mode = "D"
  Qcutoff = 30
  seqsToAnalyze = 1
#    outMode = "fq"
  numFrags = -1

  #runs the program with given inputs
  seqs = readFile(fileIn)
  seqs = rollingTrim(seqs, mode, Qcutoff, windowSize, minLen)
  fq.removeAllLeadingAndTrailingNs(seqs)
  frags = fq.seqsToFragments(seqs, minLen, numFrags)
  writeFile(fileOut, frags, mode, Qcutoff, windowSize, minLen, seqsToAnalyze)
  
  print("program end:")

