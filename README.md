# ProwlerTrimmer
Trimming tool for Oxford Nanopore sequence data

var1: [filename] = filename, including file extention (string)

var2: [inFolder] = source folder (string)

var3: [outFolder] = output folder (string)

var4: [windowSize] = number of bases in trimmer window (integer)

var5: [minLen] = minimum number of bases in read. reads with fewer bases will be rejected

var6: [trimSpecs] = trim specs: [X1]-[X2]-[X3] (string)
	X1: "L", "T", "LT", or "". L clips leading Ns, T clips trailing Ns, LT clips both, "" clips neither. recommended mode is LT
	X2: U0, F0, F1, F2... U0=unfragmented output, F0=fragmented output with all fragments, F[n]= fragmented output with n largest fragments
	X3: S, D. S=static mode. D=dynamic mode
  example: "X1-X2-X3" -> "LT-U0-S"

var7: [Qcutoff] = phred quality score threshold. (integer)

var8: [seqsToAnalyze] = megabytes of data to trim. trimmer will read files in 1 MB chunks and stop when number of MB exceeds this number.

var9: [outMode]: output file extension. saves trimmed data as either ".fasta" or ".fastq"

example execution:

generic:
python3 TrimmerLarge.py [filename] [inFolder] [outFolder] [windowSize] [minLen] [trimSpecs] [Qcutoff] [seqsToAnalyze] [outMode]

with example vriables:
python3 TrimmerLarge.py "myFile.fastq" "~/myFastqFiles/" "~/myTrimmedFiles/" 1000 1000 "LT-F1-S" 12 10 ".fastq"
