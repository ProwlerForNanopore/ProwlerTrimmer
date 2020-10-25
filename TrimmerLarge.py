# -*- coding: utf-8 -*-
"""
Created on Sun May 24 23:46:11 2020

@author: simon
"""

import Fastq as fq
import time as tm
import sys

def generateOutputFilename(filename, outFolder, trimSpecs, Qcutoff, windowSize, \
                           minLen, seqsToAnalyze, suffix):
  outputFileName = filename[0:-6] + suffix
  outputFileName += trimSpecs + str(Qcutoff)
  outputFileName += "W" 
  outputFileName += str(windowSize) 
  outputFileName += "L" + str(minLen)
  outputFileName += "R" + str(seqsToAnalyze)
  outputFileName = outFolder + outputFileName
  return outputFileName

def writeFile(filename, outFolder, seqs, trimSpecs, Qcutoff, windowSize, minLen, \
              seqsToAnalyze, suffix, outMode=".fastq"):
  # add code suffixes to output file name for mode(D__/S__), window size (W__),
  # minimum Length of read (L__), and sequences selected from raw file (R[__-__)
  
  #generates the output filename from sourcefilename and specs
  outputFileName = generateOutputFilename(filename, outFolder, trimSpecs, Qcutoff, \
                                          windowSize, minLen, seqsToAnalyze, \
                                          suffix)
  
  if outMode == ".fastq":
    fq.writeFastq(seqs,outputFileName)
  elif outMode == ".fasta":
    fq.writeFasta(seqs,outputFileName)
  else:
    print("error. output file type not supported")

# clips ends according to clipEnds seeting
def clipFunction(clipEnds,seqs):
  if clipEnds == "LT":
    fq.removeAllLeadingAndTrailingNs(seqs)
  elif clipEnds == "L":
    fq.removeAllLeadingNs(seqs)
  elif clipEnds == "T":
    fq.removeAllTrailingNs(seqs)

# Trims contig list based on min len and desired #fragments  
def fragFunction(seqs, minLen, numFrags, fragMode):
  if fragMode == "F":
    seqs = fq.seqsToFragments(seqs, minLen, numFrags)    
  return seqs

def analyzeFastqFile(fileLines,minLen):
  seqQty = int(len(fileLines)/4)
  seqs = [None]*seqQty
  for i in range(seqQty):
#    print("i:",i)
    seqs[i] = fq.Fastq(fileLines[i*4:i*4+4],minLen)
  return seqs

def readBuffer(f, bufferSize):
  fileLines = f.readlines(bufferSize)
  extraLines = (4-len(fileLines)%4)%4
  for j in range(extraLines):
    fileLines.append(f.readline())
  print(len(fileLines)/4)
    
  return fileLines

if __name__ == "__main__":
  
  #use homePC=1 when running in spyder, use 0 when running on HPC
  homePC = 1
  if homePC == 1:
    #the source file that will be trimmed
    filename = "brahSplits00-00.fastq"
    
    #folder in which the source file is located
    inFolder = ""
    
    #folder where the trimmed file will be saved
    outFolder = ""
    
    #size of the trim window
    windowSize = 1000

    #minimum length of each fastq sequence
    # minLen = 1000
    
    #trim mode (D=dynamic, S=static)
    mode = "S"
    
    #Phred Quality score threshold
    Qcutoff = 0
    
    #number of sequences to analyze (0=all)
    seqsToAnalyze = 0
    
    #u\output format of the trimmed file (fq=fastq, fa=fasta)
    outMode = ".fastq"
    
    #three trim settings:
    #[0]: leading / trailing trim flags L=leading, T=triailing
    #[1]: number of fragments to keep in fragmented mode (0=all)
    #[2]: fragmentation mode: N=no fragmentation, F=fragmented
    trimSpecs = "LT-U0-S"
  else:
    filename = sys.argv[1]
    inFolder = sys.argv[2]
    outFolder = sys.argv[3]
    windowSize = int(sys.argv[4])
    minLen = int(sys.argv[5])
    # mode = sys.argv[6][-1]
    Qcutoff = int(sys.argv[7])
    seqsToAnalyze = int(sys.argv[8])
    outMode = sys.argv[9]
    # trimStats = sys.argv[10].split(".")
    trimSpecs = sys.argv[6]
  
  trimSpecsList = trimSpecs.split("-")
  
  
  #variables denoting the size of each buffer element: n*Mb = n megabytes
  n = 4
  Mb = 1024**2
  bufferSize = n*Mb
  
  #trim mode (D=dynamic, S=static)
  mode = trimSpecsList[2]
  
  #clip ends contains the leading/trailing trim flags
  clipEnds = trimSpecsList[0]
  
  #fragMode is the fragmentation flag
  fragMode = trimSpecsList[1][0] 

  #numFrags contains the number of frags to keep in fragmented mode
  numFrags = int(trimSpecsList[1][1:])
  
  #filepath combines the source file name witht he source folder to
  #produce the source filepath
  filepath = inFolder + filename

  # performs a trim of the sequence list ofr the given trim specs
  # also removes leading and trailing Ns and breaks 
  # fragmented sequences into contigs  
  print("start:")

  #generates a new blank file with the appropriate filename
  suffix1 = "Trim"
  outputFileName = generateOutputFilename(filename, outFolder, trimSpecs, \
                                          Qcutoff, windowSize, minLen, \
                                          seqsToAnalyze, suffix1)
  f = open(outputFileName+outMode,"w+")
  f.close()  
  
  f = open(filepath, "r")
  startTime = tm.time()
  
  # reads n Mb of data, then adds between 0-3 lines to ensure the number of 
  # lines is a multiple of 4 (4 lines per fastq element)
  fileLines = readBuffer(f, bufferSize)
  
  counter=0
  time0=tm.time()
  timeCurrent=0
  maxCount = round(seqsToAnalyze/1000,0)

  #analyzes buffer seqs in lots until the End of File or max count is reached
  while (len(fileLines) > 0) and (maxCount == 0 or counter <= maxCount):
    
    counter+=1
    print("count: ", counter)
    timePrevious = timeCurrent
    timeCurrent = tm.time()-time0
    timeStep = timeCurrent-timePrevious
    timeAvg = (timeCurrent+timeStep)/(counter)
    
    print("time current: ", round(timeCurrent,2), "s")
    print("time step   : ", round(timeStep,2), "s")
    print("time average: ", round(timeAvg,2), "s")
    
    #turns line strings into fastq elements and collates them into a list
    seqs1 = analyzeFastqFile(fileLines, minLen)

    #performs rolling trim on all buffer seqs and returns coverage
    trimCoverage = fq.rollingTrimAllSeqs(seqs1,Qcutoff,windowSize,mode, minLen)
    
    # clips ends according to clipEnds setting
    clipFunction(clipEnds, seqs1)
    
    # Trims contig list based on min len and desired #fragments  
    seqs1 = fragFunction(seqs1, minLen, numFrags, fragMode)
    
    #writes current trimmed buffer seqs to file
    writeFile(filename, outFolder, seqs1, trimSpecs, Qcutoff, windowSize, minLen, \
              seqsToAnalyze, suffix1, outMode)

    # reads n Mb of data, then adds between 0-3 lines to ensure the number of 
    # lines is a multiple of 4 (4 lines per fastq element)
    fileLines = readBuffer(f, bufferSize)
    
  f.close()
  
  #calculates the time taken to analyze all seqs
  endTime = tm.time()
  trimTime = endTime - startTime
  
  #the name of the csv file that will contain the trim time and coverage stats
  fileOutName = outFolder + filename[0:-6] + "Stats"+str(trimSpecs)+str(Qcutoff) \
             +"W"+str(windowSize)+"L"+str(minLen)+"R"+str(seqsToAnalyze)+".csv"
             
  # the run specifications plus the trim time and coverage stats for the run
  trimOut = (suffix1 + "\t"+ trimSpecs +"\t"+str(Qcutoff)+"\t"+str(windowSize) \
             +"\t"+str(minLen)+"\t"+str(seqsToAnalyze)+"\t"+str(trimCoverage) \
             +"\t"+str(trimTime)+"\n")

  #outputs time statistics to csv
  fo = open(fileOutName, "w+")
  fo.write(trimOut)
  fo.close()
  print("trim time: ", round(trimTime,2), "s")
  