import sys
#to read compresed files
import gzip
import re
imputfile=input("Please, write the name of your file: ")
typefile= re.search(r'\S*\.gz', imputfile)
if typefile:
   with gzip.open(imputfile) as f:
      for line in f:
        print(line)
else:
   with open(imputfile) as f:
      for line in f:
         print(line)

#Autodetect the quality scale of a file (phred+33 or phred+64
def detect_quality():
# use python cookbook
# filter whether sanger or illumina?
# S - Sanger        Phred+33,  raw reads typically (0, 40)
# X - Solexa        Solexa+64, raw reads typically (-5, 40)
# I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
# J - Illumina 1.5+ Phred+64,  raw reads typically (3, 41)
# with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
# (Note: See discussion above).
# L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
# P - PacBio        Phred+33,  HiFi reads typically (0, 93

# Trim X nucleotides from the 3' of each read given user input
def trim3_user():

# Trim each read from the 3' based on quality, either as minimum (single residue) or mean of moving window.
def trim3_quality():

# Trim X nucleotides from the 5' of each read given user input
def trim5_user():

# These trims from 5' can be done independently of each other in above order.
def trim5_quality():

# Filter out reads with a mean quality lower than specified after trimming.
def filter_quality():

# Filter out reads that are shorter than specified after trimming.
def filter_short():

# Filter out reads that have a more than a specified number of N bases (unknown bases).
def filter_bases():
# until position 25 no errors -is our file filtered? in @ line Y=filtered,N otherwise
