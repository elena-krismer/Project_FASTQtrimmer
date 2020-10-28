import sys
#to read the fastq files
import gzip
import re

inputfile = input("Please, write the name of your file: ")
typefile = re.search(r'\S*\.gz', inputfile)
if typefile:
   with gzip.open(inputfile) as f:
      for line in f:
        print(line)
else:
   with open(inputfile) as f:
      for line in f:
         print(line)

# Maybe read file twice
# first read: triming - storing in file
# second read: filtering

# make a list with all sequences and a list with phred scores use indexing to get corresponeding mean(phred_score)


nbases_dict = dict() # creating dictionary to filter later

# main function? - perform all the functions on the lines and write into file?
# perform in which order to write into file
def first_clean(infile, trimfile):
   with open('infile', 'r') as infile:
      with open ('trimfile', 'w') as trimfile:
         for line in infile:
            if line.startswith('@'):
               outfile.write(line)
               seq_line = ''
               quality_line = ''
            if line.startswith('A', 'C', 'G', 'T', 'N'):
               trim3_user(line)
               trim5_user(line)
               #
               outfile.write(line)
               seq_line = line
            if line.startswith('+'):
               outfile.write(line)
            perform same triming on phred score
            trim3_user(line)
            trim5_user(line)
               outfile.write()
infile.close()

def second_clean(file): # +pass variables for filter functions
   with open('trimfile', 'r') as trimfile:
      with open('finalfile', 'w') as finalfile:
         for line in trimfile:
            if line.startswith('@'):
               line_one = line
            if line.startswith('A', 'C', 'G', 'T', 'N'):
               seq_line = line
               for base in seq_line:
                  if base == 'N':
                     nbases_dict[seq_line] += 1 # saving number of ns in dict
            if line.startswith('+'):
               line_three = line
            qual_line = line
            if line_one is not None and seq_line is not None and line_three is not None and qual_line is not None:
               if filter_quality(seq_line) == True and filter_bases(seq_line) == True and filter_short(seq_line, threshold_reads) == True:
                  finalfile.write(line_one + seq_line + line_three + qual_line)
            line_one, seq_line, line_three, qual_line = '', '', '', ''


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
# 5' -> 3'
# also trim phred score?
def trim3_user(seq_line, n_trim3):
   t3u_seq = seq_line[:n_trim3]
   return t3u_seq

# Trim each read from the 3' based on quality, either as minimum (single residue) or mean of moving window.
def trim3_quality():

# Trim X nucleotides from the 5' of each read given user input
def trim5_user(seq_line, n_trim5):
   t5u_seq = seq_line[n_trim5:]
   return t5u_seq

# These trims from 5' can be done independently of each other in above order.
def trim5_quality():

# Filter out reads with a mean quality lower than specified after trimming.
def filter_quality(qual_line, quality):
   # probably error because of dtype
   if mean(qual_line) > quality:
      return True
   else:
      return trim += 1
      return False
# pop sequence list when with index of phred score or make dict?

# Filter out reads that are shorter than specified after trimming.
def filter_short(seq_line, threshold_reads):
   if len(seq_line) > threshold_reads:
      return True
   else:
      return add to trimmed

# Filter out reads that have a more than a specified number of N bases (unknown bases).
def filter_bases(seq_line, n_bases):
   if nbases_dict.keys() > n_bases:
      return True
   else:
      return add to trimmed
   # use created dict and
   nbases_dict.keys().pop()
# until position 25 no errors -is our file filtered? in @ line Y=filtered,N otherwise
