# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 00:50:03 2020

@author: simon
"""
import Fastq as fq
#import pylab as pl
import matplotlib.pyplot as plt
import time as tm
import numpy as np
#import collections as cl
print("program start:")

#%%
seqsToAnalyze = 1000
minLen = 1000
windowSize = 500
#file = "brahSplits00-00TrimLT-U0-S0W1000L1000R0.fastq"
file = "brahSplits00-00.fastq"
#file = "SRR6364637_1-MesoplasmaSyrphidaeStrainYJS-MinIONRaw.fastq"


#shows the number of seqs to analyze
if seqsToAnalyze == 0 :
  seqsToAnalyze = int(fq.getFileLen(file)/4)
  print("number of seqs: ", fq.getFileLen(file)/4)
else:
  print("number of seqs: ", seqsToAnalyze)
  
#reads the fastQ file and creates a list of fastq sequences
readStart = tm.time()
seqs = fq.readFastqFile(file,seqsToAnalyze, minLen=minLen)

#%%
readQcutoff = 0
windowQcutoff = 8
#creates lists for average quality, qscore sequence, and positions
avQlist = []
qList = []
posList = []

readQlist = [0,7,9,11,13,15,17,19,21,25,29]
deltaList = [0,1,2,3,4,5,6,7,8,9,10]
readsList = np.zeros((len(deltaList),len(readQlist)),"int")


#    print("readQ:", readQ, ", delta:", delta)
print("start loop")
for n in range(seqsToAnalyze):
  seq = seqs[n]
#  print("n: ", n)
  
  #returns the rolling window quality scores for each seq
  seq.setRollingWindowQuality(windowSize)
  Q=seq.rollingWindowQuality

  #returns the positions of each rolling window start point for each seq
  P=seq.rollingWindowPos
  
  #gets the average quality for each seq
  avQ = seq.averageQuality
  Q2=sorted(Q)
  
  for ix in range(len(readQlist)):
    for iy in range(len(deltaList)):
      if len(Q2) >= 4:
        if  Q2[-2] >= (Q2[0]+deltaList[iy]) and avQ <= readQlist[ix]:
          readsList[iy,ix] += 1

    # adds avQ,P, and Q to their lists
  if len(Q2)>2:
    if avQ < windowQcutoff and Q2[-3] > windowQcutoff and len(Q) >5:
      avQlist.append(avQ)
      posList.append(P)
      qList.append(Q)
  #  
      print("Q,P,AvQ:",Q, P, avQ)
print(readsList) 
    
#%%

#fig, axs = plt.subplots(2,2)
#fig.suptitle('subtitle here')




#%%
#pl.figure(1)
#for i in range(len(avQlist)):
#  if avQlist[i] < readQcutoff:
#    avQlist[i] = 0
#axs[0,0].plot(avQlist, '+')
#
#totQ=0
#totQtrim = 0
#sumQ=0
#sumQtrim = 0
#for i in range(len(avQlist)):
#  if (avQlist[i] >= readQcutoff):
#    if (seqs[i].trimLength >= minLen):
#      totQtrim += 1
#      sumQtrim += avQlist[i]
#    if (seqs[i].length >= minLen):
#      totQ +=1
#      sumQ += avQlist[i]
#newAvg = float(sumQ)/totQ 
#newAvgTrim = float(sumQtrim/totQtrim)
#x = [0,seqsToAnalyze]
#y = [newAvg]*2
#yTrim = [newAvgTrim]*2
##pl.plot(x,y)
#axs[0,0].plot(x,y)
#axs[0,0].plot(x,yTrim)
#
#
##%%
#
##pl.figure(2)
for i in range(len(avQlist)):
  newList = [None]*len(qList[i])
  for j in range(len(qList[i])):
    if qList[i][j] != 0:
      newList[j] = qList[i][j]
#  axs[0,1].plot(posList[i],qList[i])
#  axs[1,1].plot(posList[i],newList)
#  
  ymin = [windowQcutoff]*2
  
  xmin = [0,max([max(i) for i in posList])]
#  
#  axs[0,1].plot(xmin,ymin)
#  axs[1,1].plot(xmin,ymin)
  
  plt.figure(2)
  plt.plot(posList[i],newList)
  plt.plot(xmin,ymin)
  plt.xlabel("read position")
  plt.ylabel("window average Q score")
  plt.ylim([0,15])
#  
#%% 
  
#pl.figure(3)
totalOverLen = 0
lenList = []
for i in range(len(avQlist)):
  if len(seqs[i].sequence) >= minLen:
    totalOverLen += 1
    lenList.append(len(seqs[i].sequence))
#pl.semilogy(lenList,"+")
#axs[1,0].semilogy(lenList,"+")

#%%
print("program end:")