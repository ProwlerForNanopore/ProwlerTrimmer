# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 00:50:03 2020

@author: simon
"""
import Fastq as fq
#import pylab as pl
import matplotlib.pyplot as plt
import time as tm
#import numpy as np
import collections as cl
print("program start:")

#%%
seqsToAnalyze = 400
minLen = 1000
file = "brahSplits00-00.fastq"
#file = "SRR9987849_1.fastq"

#seqsToAnalyze = int(fq.getFileLen(file)/4)
print("number of seqs: ", seqsToAnalyze)
#reads the fastQ file and creates a list of fastq sequences
readStart = tm.time()
seqs = fq.readFastqFile(file,seqsToAnalyze, minLen=minLen)
print("read time: ", tm.time()-readStart)
#    seqs = fq.readFastqFile("sequenceData.fastq,seqsToAnalyze")
#%%  
# sets up the rolling quality plot for each sequence
#for n in range(seqsToAnalyze):
#  seqs[n].setRollingWindowQuality()
#%%
# uses rolling window trimmer for each sequence
# then sets up the rolling quality plot for each sequence
windowQcutoff = 10
readQcutoff = 10
windowSize = 500
mode = "dynamic"
print("mode = ", mode)
startTime = tm.time()
fq.rollingTrimAllSeqs(seqs[0:seqsToAnalyze],windowQcutoff,windowSize,mode,minLen)
endTime = tm.time()
print("end trim. Time elapsed: ", endTime-startTime)
#%%

#contigsSingle = seqs[0].generateContigs(0)
#contigsMulti = fq.generateContigsList(seqs[0:seqsToAnalyze],1000,1)
#contigLenList = []
#for i in range(len(contigsMulti)):
#  contigLenList.append(len(contigsMulti[i]))
#print(contigLenList)
#lenListLib = cl.Counter(contigLenList)
#print(lenListLib)

#%%

fq.writeFastq(seqs)
print("read to write took: ", tm.time() - readStart)  
#fq.writeFasta(seqs)

#%%

#creates lists for average quality, qscore sequence, and positions
avQlist = [None]*seqsToAnalyze
qList = [None]*seqsToAnalyze
posList = [None]*seqsToAnalyze


for n in range(seqsToAnalyze):
  seq = seqs[n]
  
  #returns the rolling window quality scores for each seq
  Q=seq.rollingWindowQuality
  qList[n] = Q
  
  #returns the positions of each rolling window start point for each seq
  P=seq.rollingWindowPos
  posList[n] = P
  
  #gets the average quality for each seq
  avQ = seq.averageQuality

  # adds avQ to avQlist
  avQlist[n] = avQ
    
#%%
#plt.figure(1)
fig, axs = plt.subplots(2,2)
fig.suptitle('subtitle here')
#axs[0].plot(x, y)
#axs[1].plot(x, -y)
    
#%%
    
#pl.figure(1)
for i in range(len(avQlist)):
  if avQlist[i] < readQcutoff:
    avQlist[i] = 0
axs[0,0].plot(avQlist, '+')

totQ=0
totQtrim = 0
sumQ=0
sumQtrim = 0
for i in range(seqsToAnalyze):
  if (avQlist[i] >= readQcutoff):
    if (seqs[i].trimLength >= minLen):
      totQtrim += 1
      sumQtrim += avQlist[i]
    if (seqs[i].length >= minLen):
      totQ +=1
      sumQ += avQlist[i]
newAvg = float(sumQ)/totQ 
newAvgTrim = float(sumQtrim/totQtrim)
x = [0,seqsToAnalyze]
y = [newAvg]*2
yTrim = [newAvgTrim]*2
#pl.plot(x,y)
axs[0,0].plot(x,y)
axs[0,0].plot(x,yTrim)


#%%

#pl.figure(2)
for i in range(seqsToAnalyze):
  newList = [None]*len(qList[i])
  for j in range(len(qList[i])):
    if qList[i][j] != 0:
      newList[j] = qList[i][j]
  axs[0,1].plot(posList[i],qList[i])
  axs[1,1].plot(posList[i],newList)
  
  ymin = [windowQcutoff]*2
  xmin = [0,max([max(i) for i in posList])]
  
  axs[0,1].plot(xmin,ymin)
  axs[1,1].plot(xmin,ymin)
  
  plt.figure(2)
  plt.plot(posList[i],newList)
  plt.plot(xmin,ymin)
  plt.xlabel("read position")
  plt.ylabel("window average Q score")
  plt.ylim([0,30])
  
#%% 
  
#pl.figure(3)
totalOverLen = 0
lenList = [None]*seqsToAnalyze
for i in range(seqsToAnalyze):
  if len(seqs[i].sequence) >= minLen:
    totalOverLen += 1
    lenList[i] = len(seqs[i].sequence)
#pl.semilogy(lenList,"+")
axs[1,0].semilogy(lenList,"+")

#%%
print("program end:")