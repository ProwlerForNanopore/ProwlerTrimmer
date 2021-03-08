# ProwlerTrimmer
Trimming tool for Oxford Nanopore sequence data

For runMode = 0

-f, 	--file	 	filename	The name of the file you want to trim, wihtout the folderpath"

-i, 	--infolder 	inFolder	The folderpath where your file to be trimmed is located (default = cwd)

-o, 	--outfolder 	outFolder	The folderpath where your want to save the trimmed file (default = cwd)

-w, 	--windowsize 	windowSize	Change the size of the trimming window (default= 100bp)

-l, 	--minlen 	minLen		Change the minimum acceptable numer of bases in a read (default=100)

-m, 	--trimmode 	mode		Select trimming algorithm: S for static  or D for dynamic (default=S)

-q, 	--qscore 	Qcutoff		Select the phred quality score trimming threshold (default=7)

-d, 	--datamax 	maxDataMB	Select a maximum data subsample in MB (default=0, entire file)

-r, 	--outformat 	outMode		Select output format of trimmed file (fastq or fasta) (default=.fastq)

-c, 	--clip	 	clipping	Select L to clip leading Ns, T to trim trialing Ns and LT to trim both (default=LT)

-g, 	--fragments 	fragments	Select fragmentation mode (default=U0)


[filename] = filename, including file extention (string)

[inFolder] = source folder (string)

[outFolder] = output folder (string)

[windowSize] = number of bases in trimmer window (integer)

[minLen] = minimum number of bases in read. reads with fewer bases will be rejected

[trimSpecs] = trim specs: [X1]-[X2]-[X3] (string)
	X1: "L", "T", "LT", or "". L clips leading Ns, T clips trailing Ns, LT clips both, "" clips neither. recommended mode is LT
	X2: U0, F0, F1, F2... U0=unfragmented output, F0=fragmented output with all fragments, F[n]= fragmented output with n largest fragments
	X3: S, D. S=static mode. D=dynamic mode
  example: "X1-X2-X3" -> "LT-U0-S"

[Qcutoff] = phred quality score threshold. (integer)

[seqsToAnalyze] = megabytes of data to trim. trimmer will read files in 1 MB chunks and stop when number of MB exceeds this number.

[outMode]: output file extension. saves trimmed data as either ".fasta" or ".fastq"

example execution:



**runMode=0:**
python3 TrimmerLarge.py -f [filename] -i [inFolder] -o [outFolder] -w [windowSize] -l [minLen] -c [clipping] -g [fragments] -m [mode] -q [Qcutoff] -d [maxDataMB] -r [outMode]

with example vriables:
python3 TrimmerLarge.py -f "myFile.fastq" -i "~/myFastqFiles/" -o "~/myTrimmedFiles/" -w 1000 -l 1000 -c "LT" -g "F1" -m "S" -q 12 -d 10 -r ".fastq"

or with defaults, selecting only the algorithm mode and the quality cutoff:
python3 TrimmerLarge.py -f "myFile.fastq" -m "S" -q 12

**runMode=2:**
This mode requires the variables be listed in the given order and recalled using sys.arg[n]. 
This usage is deprecated and currently runMode is a hardcoded variable that must be changed by editing the python file. 
The preferred way to run this script is using runMode=0.

python3 TrimmerLarge.py [filename] [inFolder] [outFolder] [windowSize] [minLen] [trimSpecs] [Qcutoff] [seqsToAnalyze] [outMode]

with example vriables:
python3 TrimmerLarge.py "myFile.fastq" "~/myFastqFiles/" "~/myTrimmedFiles/" 1000 1000 "LT-F1-S" 12 10 ".fastq"

**runMode=1:**
A mode that requires the variables section of the python script to be set by editing the python file.

