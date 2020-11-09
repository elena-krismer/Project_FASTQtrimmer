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

# global variables

# user input
nt3_input, nt5_input, threshold_reads_input, quality_input, n_bases_input = None, None, None, None, None

# counting triming

# filtered reads is this possible as global variable????
filter_quality_count: int
filter_quality_count, filter_short_count, filter_bases_count = 0, 0, 0


# make a list with all sequences and a list with phred scores use indexing to get corresponeding mean(phred_score)


nbases_dict = dict() # creating dictionary to filter later
#Autodetect the quality scale of a file (phred+33 or phred+64)
def detect_quality(ascii_string):
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
      return filter_quality_count += 1

# pop sequence list when with index of phred score or make dict?

# Filter out reads that are shorter than specified after trimming.
def filter_short(seq_line, threshold_reads):
   if len(seq_line) > threshold_reads:
      return True
   else:
      return filter_short_count += 1

# Filter out reads that have a more than a specified number of N bases (unknown bases).
def filter_bases(seq_line, n_bases):
   if nbases_dict.keys() > n_bases:
      return True
   else:
      return filter_bases_count += 1 # keeping track


# main function? - perform all the functions on the lines and write into file?
# perform in which order to write into file
def first_clean(infile, trimfile):
   with open('infile', 'r') as infile:
      with open ('trimfile', 'w') as trimfile:
         for line in infile:
            if line.startswith('@'):
               trimfile.write(line)
            elif re.search(r^'(ACGTN)', line) is None:
               # make trimint function???
               line = trim3_user(line, nt3_input)
               line = trim5_user(line, nt3_input)
               line = trim3_quality(line, XX)
               line = trim5_quality(line, XX)
               trimfile.write(line)
            elif line.startswith('+'):
               trimfile.write(line)
            elif: # perform same trimming on phred score
               phred_scale = detect_quality(line)
               line = trim3_user(line, nt3_input)
               line = trim5_user(line, nt5_input)
               line = trim3_quality(line, XX)
               line = trim5_quality(line, XX)
               trimfile.write(line)
            else:
               print('Error in Fileformat.')
               sys.exit(1)
infile.close()

def second_clean(trimfile, finalfile): # +pass variables for filter functions
   with open('trimfile', 'r') as trimfile:
      with open('finalfile', 'w') as finalfile:
         line_one, seq_line, line_three, qual_line = '', '', '', ''
         for line in trimfile:
            if line.startswith('@'):
               line_one = line
            if re.search(r'^(ACGTN)', line) is None:
                  seq_line = line
                  for base in seq_line:
                     if base == 'N' :
                        nbases_dict[seq_line] += 1 # saving number of ns in dict
            if line.startswith('+'):
                  line_three = line
                  qual_line = line +1

            # check if all filters are passed than write into file and set lines to None
            if line_one is not None and seq_line is not None and line_three is not None and \
               qual_line is not None:
               if filter_quality(qual_line, quality_input) == True and \
                  filter_bases(nbases_dict[seq_line], n_bases_input) == True and \
                  filter_short(seq_line, threshold_reads_input) == True:
                     finalfile.write(line_one + '\n' + seq_line + '\n' + line_three + '\n' + qual_line)
            line_one, seq_line, line_three, qual_line = '', '', '', ''

def summary_file(summaryfile):
   with open ('summaryfile', 'w') as sum_file:
      filtered = filter_quality_count + filter_bases_count + filter_short_count
      sum_file.write('Total number of reads filtered: ' + filtered + filter_quality_count + \
                     'reads with a low quality than' + quality_input + filter_bases_count + 'reads with more than' \
                     + filter_bases_count + 'bases')

# use created dict and
   nbases_dict.keys().pop()
# until position 25 no errors -is our file filtered? in @ line Y=filtered,N otherwise
