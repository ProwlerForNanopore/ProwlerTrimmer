# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 00:26:41 2020

@author: simon
"""
import numpy as np
#import pylab as pl

class Fastq:
  #initializes the class by taking the four lines of the fastQ format
  #and turns that into (0)name, (1)sequence, and (3)quality
  def __init__(self, fourLines,minLen):
    self.name = fourLines[0]
    self.qualityASCII = fourLines[3].strip()
    self.length = len(self.qualityASCII)
    
    if self.length >= minLen:
      self.discard = 0
      self.trimLength = self.length
      self.sequenceStr = fourLines[1].strip()
      self.sequence = [None]*self.length
      for i in range(self.length):
        self.sequence[i] = self.sequenceStr[i]
  
      self.quality = np.zeros(self.length, "int") 
      for i in range(self.length): 
        self.quality[i] = ord((self.qualityASCII[i]))-33
      self.averageQuality = sum(self.quality)/self.length
  
      self.rollingWindowQuality = np.array([], "int") #<--------------------------------------------modified: None
      self.rollingWindowPos = np.array([], "int") #<--------------------------------------------modified: None
      self.windowSize = 0 #<--------------------------------------------modified: None
      self.fastqList = []
      self.minLen = minLen
    else:
      self.discard = 1
      self.length = 0
      self.trimLength = 0
      self.sequenceStr = ""
      self.sequence = []
      self.quality = np.array([])
      self.averageQuality = 0
      self.rollingWindowQuality = np.array([], "int") #<--------------------------------------------modified: None
      self.rollingWindowPos = np.array([], "int") #<--------------------------------------------modified: None
      self.windowSize = 0 #<--------------------------------------------modified: None
      self.fastqList = []
      self.minLen = minLen
    
    
  #for the purpose of sorting, a less than comparison function is required
  def __lt__(self, other):
    return self.trimLength < other.trimLength
  
  #After trimming, all trimmed sequences are set to N and all trimmed qualitys 
  #are set to 0. When converting back to ASCII, this method redefines the ASCII 
  #string before outputting to file.
  def setQualityASCII(self):
    if self.discard == 1: return
    temp = [None]*len(self.quality)
    for i in range(len(self.quality)):
      temp[i] = chr((self.quality[i])+33)
    self.qualityASCII = "".join(temp)
  
  
  #After trimming, the average quality is recalculated by taking the score of 
  #all remaining elements and dividing it by the number of remaining elements.
  def setAverageQuality(self, windowSize):
    if self.discard == 1: return
    sumQ=0
    totQ=0
    numberOfWindows = len(self.rollingWindowQuality)
    overhangLength = windowSize*numberOfWindows - self.length
    
    for i in range(numberOfWindows-1):
      #ads the score and count for each nonzero window to the quality average
      Qscore = self.rollingWindowQuality[i]
      if Qscore != 0:
        sumQ += Qscore
        totQ += 1
        
    #ads the overhand score and fractional count to the quality average
    Qscore = self.rollingWindowQuality[-1]
    if Qscore != 0:
      fractionalCount = overhangLength/windowSize
      sumQ += Qscore*fractionalCount
      totQ += fractionalCount
    
    #calculates the average quality only if the window count is nonzero
    if totQ == 0:
      self.averageQuality = 0
    else :
      self.averageQuality = float(sumQ) / float(totQ)

  #Before trimming using a rolling window, the rolling windows must be defined. 
  #This method creates a list of window qualities and positions. Window 
  #quality, position and size are given values in this method.
  def setRollingWindowQuality(self,windowSize):
    if self.discard == 1: return
    
    self.windowSize = windowSize
    
    #calculates number of windows
    numberOfWindows = int(self.length/windowSize + 1)
    
    #creates lists for each window position and score
    windowsPos = np.zeros(numberOfWindows,"int")
    scores = np.zeros(numberOfWindows,"int")
    
    # calculates window quality
    for windowNum in range(numberOfWindows):
      pos = windowNum*windowSize
      if pos+windowSize <= self.length:
        windowAvg = sum(self.quality[pos : pos + windowSize]) / windowSize
      elif self.length-pos > 0:
        windowAvg = sum(self.quality[pos : self.length]) / (self.length-pos)
      scores[windowNum] = windowAvg
      windowsPos[windowNum] = pos
      
    #updates class variables for rolling window information
    self.rollingWindowQuality = scores
    self.rollingWindowPos = windowsPos

  #checks if fastq element has a sequence greater than minLen
  #sequence is set to be discarded if trim length < min length
  def minLengthCheck(self):
    if self.discard == 1: return
    if self.trimLength < self.minLen:
      self.discard = 1

  #This method takes the rolling window variables and performs a trim on any 
  #window that falls below the minQ cut-off. Average quality is updated after 
  #the trim.
  def rollingTrim(self,minQ):
    if self.discard == 1: return
    
    windowSize = self.windowSize

    #calculates number of windows
    numberOfWindows = len(self.rollingWindowPos)
    
    # calculates avg quality of each window and trims all below cutoff
    for windowNum in range(numberOfWindows):
      
      #performs min length check on sequence. ends trimmer if below min len
      self.minLengthCheck()
      if self.discard == 1:
        break
      
      pos = self.rollingWindowPos[windowNum]
      if pos+windowSize <= self.length:
        windowAvg = self.rollingWindowQuality[windowNum]
        
        # trims the window if the avg q score is below cutoff
        if windowAvg < minQ and windowAvg != 0:
          for base in range(pos,pos+windowSize):
            self.sequence[base] = "N"
            self.quality[base] = 0
          windowAvg = 0
          # reduces the effective length by the trimmed window size
          self.trimLength -= windowSize

      elif self.length-pos > 0:
        windowAvg = self.rollingWindowQuality[windowNum]
        if windowAvg < minQ and windowAvg != 0:
          for base in range(pos,self.length):
            self.sequence[base] = "N"
            self.quality[base] = 0
          windowAvg = 0
          # reduces the effective length by the trimmed window size
          self.trimLength -= self.length-pos

      self.rollingWindowQuality[windowNum] = windowAvg

    # updates class instance: window quality and poistion information
    self.setAverageQuality(windowSize)

  #Static trim executed the rolling trim method a single time using minQ as 
  #the window-wide score cut-off.
  def rollingTrimStatic(self, minQ, windowSize):
    if self.discard == 1: return
    self.setRollingWindowQuality(windowSize)
    self.rollingTrim(minQ)
    
  #Dynamic trim executes the rolling window trim multiple times on each read 
  #with successively higher window-wide cut-offs until a given average read 
  #quality is obtained. The method starts with a window-wide cut-off of 2 and 
  #increases by 2 each time until either the read average is accepted or the 
  #read length falls below minReadLength (TO DO).
  def rollingTrimDynamic(self,readQmin,windowSize):
    if self.discard == 1: return
    self.setRollingWindowQuality(windowSize)
    minQ = max([readQmin-10, 2])
    self.rollingTrim(minQ)
    while (self.averageQuality < readQmin) and (self.trimLength > self.minLen):
      minQ += 2
      self.rollingTrim(minQ)

  #a helper method that adds fragments to a list only if they are above minLen
  def addContigToList(self, tempFastqList,contigName, tempSeq, tempScore):
    if len(tempSeq) >= self.minLen:
      tempFastqList.append(FastqFromVars(contigName, tempSeq, tempScore, self.minLen))

  #This method takes a trimmed, gappy fastq object and breaks it into multiple 
  #smaller fastq objects with no gaps. It rejects any sequence that is shorter 
  #than minLen. Finally, it ranks sequences by length and only accepts a given 
  #number of them for output as a list of fastq objects.
  def generateContigs(self, numberOfContigs):
    if self.discard == 1: return
    
    #sets up contig lists
    tempSeq = []
    tempScore = ""
    fastqList = []
    tempFastqList = []
#    contigName = ""
    
    #counters that will track if we are in a contig or a trim
    contigNumber = -1
    trimCount = 1
    
    #calculates number of windows
    numberOfWindows = len(self.rollingWindowPos)
    
    #main loop
    for i in range(numberOfWindows):
      
      #adds window to current/new contig depending on value of trimCount
      if self.rollingWindowQuality[i] != 0:
        
        #determines if this is a new or existing contig
        if trimCount !=0:
          tempSeq = []
          tempScore = ""
          contigNumber += 1
          trimCount = 0
          contigName = self.name[:-1] + " @window " + str(i) + "/" + \
                                            str(numberOfWindows) + "\n"
        #window size is 'self.windowSize' except wiht final window          
        start = self.rollingWindowPos[i]
        if i+1 == numberOfWindows:
          end = len(self.sequence)
        else:
          end = start + self.windowSize
          
        # updates the conig sequence and score
        tempSeq += self.sequence[start:end]
        tempScore += self.qualityASCII[start:end]
        
        #writes to fastq if this is the last window
        if i+1 == numberOfWindows:
          self.addContigToList(tempFastqList,contigName, tempSeq, tempScore)
#          tempFastqList.append(FastqFromVars(contigName, tempSeq, tempScore, self.minLen))
      
      # toggles trim region identifyer
      elif trimCount == 0:
        trimCount += 1
        self.addContigToList(tempFastqList,contigName, tempSeq, tempScore)
#        tempFastqList.append(FastqFromVars(contigName, tempSeq, tempScore, self.minLen))
          
    #removes any fragments shorter than minLen
    for fragment in tempFastqList:
      if fragment.trimLength >= self.minLen:
        fastqList.append(fragment)
    
    #sorts the list of fragments by length, 
    #then, samples qty of ragments of up to numberOfContigs
    fastqList.sort(reverse=True)
    if numberOfContigs > 0:
      if len(fastqList) > numberOfContigs:
        fastqList = fastqList[:numberOfContigs]
        
    self.fastqList = fastqList
    return fastqList

  # removes leading Ns from a fastq element. Updates RollingWindowPos, 
  # rollingWindowQuality, sequence, quality, and ASCIIquality
  def removeLeadingNs(self):
    if self.discard == 1: return
    
    # counts the number of windows on the leading edge that have a window 
    # quality of 0 (which means it has been trimmed)
    deleteLWindowCount = 0
    for i in range(len(self.rollingWindowQuality)):
      if self.rollingWindowQuality[i] == 0:
        deleteLWindowCount += 1
      else:
        break
    
    # uses deleteLWindowCount to remove all leading Q=0 windows, then deletes
    # sequences and qualities in the index range of those windows
    # Finally, reduces the value of all elements in self.rollingWindowPos
    # so that the window pos values start from zero.
    self.rollingWindowQuality = self.rollingWindowQuality[deleteLWindowCount:]
    self.rollingWindowPos = self.rollingWindowPos[deleteLWindowCount:]
    self.rollingWindowPos -= deleteLWindowCount*self.windowSize
    self.sequence = self.sequence[deleteLWindowCount*self.windowSize:]
    self.quality = self.quality[deleteLWindowCount*self.windowSize:]
    
  # removes trailing Ns from a fastq element. Updates RollingWindowPos, 
  # rollingWindowQuality, sequence, quality, and ASCIIquality
  def removeTrailingNs(self):
    if self.discard == 1: return
    
    #counts the number of windows on the trailing end that have a quality of 0
    deleteTWindowCount = 0
    for i in range(len(self.rollingWindowQuality)):
      j = len(self.rollingWindowQuality) - 1 - i
      if self.rollingWindowQuality[j] == 0:
        deleteTWindowCount += 1
      else:
        break

    # uses deleteTWindowCount to remove all trailing Q=0 windows, then dletes
    # sequences and qualities in the index range of those windows.
    if deleteTWindowCount > 0:
      #must factor in partial window on trailing edge
      overhang = self.length%self.windowSize
      deleteQty = (deleteTWindowCount-1)*self.windowSize + overhang
      
      self.rollingWindowQuality = self.rollingWindowQuality[:-deleteTWindowCount]
      self.rollingWindowPos = self.rollingWindowPos[:-deleteTWindowCount]
      self.sequence = self.sequence[:-deleteQty]
      self.quality = self.quality[:-deleteQty]
        
#------------------------------------------------------------------------------

# takes as input a list of fastq sequence elements containing N gaps.
# Breakes each individual sequence into contigs, then adds all contigs 
# from all sequences to a new, expanded list of gapless fragments.
def seqsToFragments(seqs, minLen=0, numContigs=-1):
  fragmentsList = []
  for seq in seqs:
    if seq.discard == 1: continue
    fragments = seq.generateContigs(numContigs)
    for fragment in fragments:
      fragmentsList.append(fragment)
      
  return fragmentsList
    

# removes leading Ns from a list of fastq elements.
def removeAllLeadingNs(seqs):
  for i in range(len(seqs)):
    seqs[i].removeLeadingNs()

# removes leading Ns from a list of fastq elements.
def removeAllTrailingNs(seqs):
  for i in range(len(seqs)):
    seqs[i].removeTrailingNs()
    
#removes all leading and trailing Ns from a list of fasq elements
def removeAllLeadingAndTrailingNs(seqs):
  removeAllLeadingNs(seqs)
  removeAllTrailingNs(seqs)

#This method performs a rolling trim for a given window size on every fastq 
#in a list of seqs. The mode can be either static or dynamic. Qcutoff will be 
#read-wide or window-wide depending on mode.
def rollingTrimAllSeqs(seqs,Qcutoff,windowSize,mode, minLen):
  trimCoverage = 0
  rawCoverage = 0
  if mode == "S":
#    print("entering static mode:")
    for n in range(len(seqs)):
      seqs[n].rollingTrimStatic(Qcutoff,windowSize)
      trimCoverage += seqs[n].trimLength
      rawCoverage += seqs[n].length
      
  elif mode == "D":
#    print("entering dynamic mode:")
    for n in range(len(seqs)):
      seqs[n].rollingTrimDynamic(Qcutoff,windowSize)
      trimCoverage += seqs[n].trimLength
      rawCoverage += seqs[n].length
  else:
    print("error. mode selected not available. please select S or D")
  
#  if rawCoverage > 0:
#    print("raw: ", rawCoverage, " Trim: ", trimCoverage, \
#        "ratio: ",round(trimCoverage/rawCoverage,2))
#  else:
#    print("raw coverage: ", rawCoverage)
#    print("trim coverage: ", trimCoverage)
  
  return trimCoverage

#This method generates a fastq class element from three strings representing 
#name, sequence, and Q-score
def FastqFromVars(name,seq,score,minLen):
  seqString = "".join(seq) + " "
  fourLines = [name,seqString,"+",score]
  return Fastq(fourLines, minLen) 
  
#This method executed the generate contigs method for every fastq element in 
#a list of seqs. The resulting output will be a list of fastq lists.
def generateContigsList(seqs, minLen, numberOfContigs):
  contigsList = [None]*len(seqs)
  for i in range(len(seqs)):
    contigsList[i] = (seqs[i].generateContigs(numberOfContigs))
  
  return contigsList
  
#Takes a fastq file and converts it to a list of fastq objects as defined 
#by the fastq class
def readFastqFile(FQfile,numberOfSeqs=0, coverage=0, minLen=1000):
  f = open(FQfile, "r")
  file = f.readlines()
  f.close()
  
  #if number of seqs is not given, read all seqs
  if numberOfSeqs <= 0:
    seqQty = int(len(file)/4)
  else:
    seqQty = numberOfSeqs
    
  #if coverage is given, read until coverage is reached
  if coverage > 0 : 
    seqs = []
    totalCoverage = 0
    i = 0
    while totalCoverage < coverage:
      seqs.append(Fastq(file[i*4:i*4+4]),minLen)
      totalCoverage += seqs[i].length
      i+=1
  else:
    seqs = [None]*seqQty
    for i in range(seqQty):
      seqs[i] = Fastq(file[i*4:i*4+4], minLen)
  return seqs


#Takes a list of fastq objects as defined by the fastq class and outputs a 
#fastq file in fastq fromat
def writeFastq(seqs,file):
  dcount=0
  acount=0
  file += ".fastq"
  f = open(file,"a+")
  for i in range(len(seqs)):
    if seqs[i].discard == 0:
      f.write(seqs[i].name)
      # writes sequence. takes list of strings and makes string, then writes
      seq = "".join(seqs[i].sequence)
      f.write(seq)
      
      #new liine, then a + between seq and quality, then new line again
      f.write("\n+\n")
      
      #then quality, then finally, a new line
      seqs[i].setQualityASCII()
      f.write(str(seqs[i].qualityASCII))
      f.write("\n")
      acount+=1
    else:
      dcount+=1
      
  # this if statement is to prevent div by zero error when printing the discard ratio
#  if (dcount + acount == 0):
#    print("dcount + acount = 0")
#  else:
#    print("discarded:", round(dcount/(dcount+acount),2) , \
#          "printed:", round(acount/(dcount+acount),2))
  f.close()
  
#Takes a list of fastq objects as defined by the fastq class and outputs a 
#fasta file in fasta fromat
def writeFasta(seqs, file):
  dcount=0
  acount=0
  file += ".fasta"
  f = open(file,"a+")
  for i in range(len(seqs)):
    if seqs[i].discard == 0:
      f.write(">")
      f.write(seqs[i].name)
      # writes sequence. takes list of strings and amkes string, then writes
      seq = "".join(seqs[i].sequence)
      f.write(seq)
      f.write("\n")
      acount+=1
    else:
      dcount+=1
#  print("discarded:", dcount ,"printed:", acount)
  f.close()
    
def getFileLen(FQfile):
  f = open(FQfile, "r")
  file = f.readlines()
  f.close()
  return len(file)
    
if __name__ == "__main__":
  seqs = readFastqFile("G:/Honours Thesis Data/SequenceData/testSeq.fastq", \
                       numberOfSeqs=1000)
  minLen = 1000
  mode = "S"
  windowSize = 500
  Qcutoff = 20
  trimCoverage = rollingTrimAllSeqs(seqs,Qcutoff,windowSize,mode, minLen)
  seqs = readFastqFile("SRR6364637_1-MesoplasmaSyrphidaeStrainYJS-MinIONRaw.fastq", \
                       coverage=trimCoverage)
  trimCoverage0 = rollingTrimAllSeqs(seqs,0,windowSize,mode, minLen)  
