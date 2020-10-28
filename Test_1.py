import sys
#to read the fastq files
import gzip
import re
inputfile=input("Please, write the name of your file: ")
typefile= re.search(r'\S*\.gz', inputfile)
if typefile:
   with gzip.open(inputfile) as f:
      for line in f:
        print(line)
else:
   with open(inputfile) as f:
      for line in f:
         print(line)

#Autodetect the quality scale of a file (phred+33 or phred+64
def detect_quality(ascii_string):
   phred_scale = ''
   ascii_list = [ord(ascii_value) for ascii_value in ascii_string]
   # phred 33 range= 33-75
   if max(ascii_list) <= 75 and min(ascii_list) < 59:
      phred_scale = 33
   # phred 64 range = 64 - 106
   elif max(ascii_list) > 75 and min(ascii_list) >= 64:
      phred_scale = 64
   else:
      print("Error in determining quality scale")
      sys.exit(1)
   return phred_scale

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
