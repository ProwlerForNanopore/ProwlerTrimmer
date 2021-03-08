# -*- coding: utf-8 -*-
"""
Created on Sun May 24 23:46:11 2020

@author: simon
"""

import Fastq as fq
import time as tm
import sys
from optparse import OptionParser

def generateOutputFilename(filename, outFolder, trimSpecs, Qcutoff, windowSize, \
                           minLen, maxDataMB, suffix):
  outputFileName = filename[0:-6] + suffix
  outputFileName += trimSpecs + str(Qcutoff)
  outputFileName += "W" 
  outputFileName += str(windowSize) 
  outputFileName += "L" + str(minLen)
  outputFileName += "R" + str(maxDataMB)
  outputFileName = outFolder + outputFileName
  return outputFileName

def writeFile(filename, outFolder, seqs, trimSpecs, Qcutoff, windowSize, minLen, \
              maxDataMB, suffix, outMode=".fastq"):
  # add code suffixes to output file name for mode(D__/S__), window size (W__),
  # minimum Length of read (L__), and sequences selected from raw file (R[__-__)
  
  #generates the output filename from sourcefilename and specs
  outputFileName = generateOutputFilename(filename, outFolder, trimSpecs, Qcutoff, \
                                          windowSize, minLen, maxDataMB, \
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
#  print(len(fileLines)/4)
    
  return fileLines

if __name__ == "__main__":
  
  #use homePC=1 when running in spyder, use 0 when running on HPC
  runMode = 0
  
  if runMode == 0:
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                  help="the name of the file you want to trim, wihtout the folderpath", metavar="FILE")
    parser.add_option("-i", "--infolder", dest="inFolder",
                  help="the folderpath where your file to be trimmed is located (default = cwd)", metavar="DIRECTORY", default="")
    parser.add_option("-o", "--outfolder", dest="outFolder",
                  help="the folderpath where your want to save the trimmed file (default = cwd)", metavar="DIRECTORY", default="")
    parser.add_option("-w", "--windowsize", dest="windowSize",
                  help="change the size of the trimming window (default= 100bp)", metavar="INT", default=100)
    parser.add_option("-l", "--minlen", dest="minLen",
                  help="change the minimum acceptable numer of bases in a read (default=100)", metavar="INT", default=100)
    parser.add_option("-m", "--trimmode", dest="mode",
                  help="select trimming algorithm: S for static  or D for dynamic (default=)", metavar="[S/D]", default="S")
    parser.add_option("-q", "--qscore", dest="Qcutoff",
                  help="select the phred quality score trimming threshold (default=7)", metavar="INT", default=7)
    parser.add_option("-d", "--datamax", dest="maxDataMB",
                  help="select a maximum data subsample in MB (default=0, entire file))", metavar="INT", default=0)
    parser.add_option("-r", "--outformat", dest="outMode",
                  help="select output format of trimmed file (fastq or fasta) (default=.fastq)", metavar="[.fasta/.fastq]", default=".fastq")
    parser.add_option("-c", "--clip", dest="clipping",
                  help="select L to clip leading Ns, T to trim trialing Ns and LT to trim both (default=LT)", metavar="[L/T/LT]", default="LT")
    parser.add_option("-g", "--fragments", dest="fragments",
                  help="select fragmentation mode (default=U0)", metavar="[U0/F1/F2...]", default="U0")
    (options, args) = parser.parse_args()
    
    #the source file that will be trimmed
    filename = options.filename
    
    #folder in which the source file is located
    inFolder = options.inFolder
    
    #folder where the trimmed file will be saved
    outFolder = options.outFolder
    
    #size of the trim window
    windowSize = int(options.windowSize)

    #minimum length of each fastq sequence
    minLen = int(options.minLen)
    
    #trim mode (D=dynamic, S=static)
    mode = options.mode
    
    #Phred Quality score threshold
    Qcutoff = int(options.Qcutoff)
    
    #number of sequences to analyze (0=all)
    maxDataMB = int(options.maxDataMB)
    
    #u\output format of the trimmed file (fq=fastq, fa=fasta)
    outMode = options.outMode
    
    trimSpecs = options.clipping + "-" + options.fragments + "-" + options.mode
    trimSpecsList = [options.clipping, options.fragments,options.mode]
    
  elif runMode == 1:
    #the source file that will be trimmed
    filename = "brahSplits00-00.fastq"
    
    #folder in which the source file is located
    inFolder = ""
    
    #folder where the trimmed file will be saved
    outFolder = ""
    
    #size of the trim window
    windowSize = 1000

    #minimum length of each fastq sequence
    minLen = 1000
    
    #trim mode (D=dynamic, S=static)
    mode = "S"
    
    #Phred Quality score threshold
    Qcutoff = 0
    
    #number of sequences to analyze (0=all)
    maxDataMB = 100
    
    #u\output format of the trimmed file (fq=fastq, fa=fasta)
    outMode = ".fastq"
    
    #three trim settings:
    #[0]: leading / trailing trim flags L=leading, T=triailing
    #[1]: number of fragments to keep in fragmented mode (0=all)
    #[2]: fragmentation mode: N=no fragmentation, F=fragmented
    trimSpecs = "LT-U0-S"
    trimSpecsList = trimSpecs.split("-")
  elif runMode == 2:
    filename = sys.argv[1]
    inFolder = sys.argv[2]
    outFolder = sys.argv[3]
    windowSize = int(sys.argv[4])
    minLen = int(sys.argv[5])
    # mode = sys.argv[6][-1]
    Qcutoff = int(sys.argv[7])
    maxDataMB = int(sys.argv[8])
    outMode = sys.argv[9]
    # trimStats = sys.argv[10].split(".")
    trimSpecs = sys.argv[6]
    trimSpecsList = trimSpecs.split("-")
  
  
  
  #variables denoting the size of each buffer element: n*Mb = n megabytes
  n = 1
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
                                          maxDataMB, suffix1)
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
  maxCount = round(maxDataMB,0)

  #analyzes buffer seqs in lots until the End of File or max count is reached
  while (len(fileLines) > 0) and (maxCount == 0 or counter <= maxCount):
    print(counter)
    counter+=1
#    print("count: ", counter)
    timePrevious = timeCurrent
    timeCurrent = tm.time()-time0
    timeStep = timeCurrent-timePrevious
    timeAvg = (timeCurrent+timeStep)/(counter)
    
#    print("time current: ", round(timeCurrent,2), "s")
#    print("time step   : ", round(timeStep,2), "s")
#    print("time average: ", round(timeAvg,2), "s")
    
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
              maxDataMB, suffix1, outMode)

    # reads n Mb of data, then adds between 0-3 lines to ensure the number of 
    # lines is a multiple of 4 (4 lines per fastq element)
    fileLines = readBuffer(f, bufferSize)
    
  f.close()
  
  #calculates the time taken to analyze all seqs
  endTime = tm.time()
  trimTime = endTime - startTime
  
  #the name of the csv file that will contain the trim time and coverage stats
  fileOutName = outFolder + filename[0:-6] + "Stats"+str(trimSpecs)+str(Qcutoff) \
             +"W"+str(windowSize)+"L"+str(minLen)+"R"+str(maxDataMB)+".csv"
             
  # the run specifications plus the trim time and coverage stats for the run
  trimOut = (suffix1 + "\t"+ trimSpecs +"\t"+str(Qcutoff)+"\t"+str(windowSize) \
             +"\t"+str(minLen)+"\t"+str(maxDataMB)+"\t"+str(trimCoverage) \
             +"\t"+str(trimTime)+"\n")

  #outputs time statistics to csv
  fo = open(fileOutName, "w+")
  fo.write(trimOut)
  fo.close()
#  print("trim time: ", round(trimTime,2), "s")
  
  
