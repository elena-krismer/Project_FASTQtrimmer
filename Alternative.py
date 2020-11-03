import sys
import statistics
import gzip
import re
import argparse

#I think the first step would be open the file. I did some changes because the code didn't run but now it works. In the modifications you added I saw that you suggest to save the file in a list, it's added

parser = argparse.ArgumentParser()
parser.add_argument('file')
args = parser.parse_args()

def open_file(filename):
   file_list=list()
   file_listb=list()
   typefile = re.search(r'\S*\.gz', args.file)
   if typefile:
      with gzip.open(args.file) as f:
         for line in f:
          file_list.append(line) #is the best way to save the lines??
      print(file_list) #we can remove it I just printed to see if it works
   else:
      with open(args.file) as f:
         for line in f:
          file_listb.append(line)
      print(file_listb) #we can remove it I just printed to see if it works

fun=open_file(args.file) #lines added to run the example
fun #line added to run the example




# user input - global or not  - bad for runtime
# nt3_input, nt5_input, threshold_reads_input, quality_input, n_bases_input = None, None, None, None, None

# counting triming

# filtered reads is this possible as global variable????
from pip._vendor.certifi.__main__ import args

filter_quality_count, filter_short_count, filter_bases_count = 0, 0, 0


def trim_user(seq_line, trim3, trim5):
    trim_line = seq_line[:-trim3]
    trim_line = trim_line[trim5:]
    return trim_line

# adding function for trimming window???
# problems in determining phred scale  - use maybe elif to check the next read for quality if the first doesnt work
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


def filter_quality(qual_line, quality):
    # converting to ascii values
    ascii_list = [ord(ascii_value) for ascii_value in qual_line]
    if statistics.mean(ascii_list) > quality:
        return True
    else:
        global filter_quality_count
        filter_quality_count += 1
        return filter_quality_count


def filter_bases(seq_line, n_bases):
    n_number = 0
    for base in seq_line:
        if base == 'N':
            n_number += 1
    if n_number < n_bases:
        return True
    else:
        global filter_bases_count
        filter_bases_count += 1
        return filter_bases_count


def filter_short(seq_line, threshold_reads):
    if len(seq_line) > threshold_reads:
        return True
    else:
        global filter_short_count
        filter_short_count += 1
        return filter_short_count


def summary_file(summaryfile, quality, nbases, length):
    with open('summaryfile', 'w') as sum_file:
        global filter_quality_count, filter_bases_count, filter_short_count
        filtered = filter_quality_count + filter_bases_count + filter_short_count
        sum_file.write(
            "Total number of reads filtered: {0}. {1} reads with a low quality than {2}. {3} reads with more \
            than {4} bases. {5} reads shorter than {6} nucleotides.".format(filtered,
                                                                            filter_quality_count, quality,
                                                                            filter_bases_count, nbases,
                                                                            filter_short_count, length))


def run(args):
    file_list = list()
    # nputfile = input("Please, write the name of your file: ")
    # typefile = re.search(r'\S*\.gz', inputfile)
    # if typefile:
    #     with gzip.open(inputfile) as f:
    #         for line in f:
    #             print(line)
    # else:
    #     with open(inputfile) as f:
    #         for line in f:
    #             print(line)
    # typefile = re.search(r'\S*\.gz', args.input)  # works filedirectory as string?
    # if typefile:
    #     infile = gzip.open(args.input)
    # else:
    #     infile = open(args.input)
    with open(args.input, 'r') as infile:
        for line in infile:
            file_list.append(line)
        infile.close()

    # trimming sequence
    pos_seq = 1
    while pos_seq < len(file_list):
        file_list[pos_seq] = trim_user(file_list[pos_seq], args.trim3, args.trim5)
        pos_seq += 4

    # trimming quality scale
    pos_qual = 3
    while pos_qual < (len(file_list) - 1):  # -1 or not?
        file_list[pos_qual] = trim_user(file_list[pos_qual], args.trim3, args.trim5)
        pos_qual += 4
    q_scale = detect_quality(file_list[6])

    pos = 0
    outputfile = open(args.output)
    while pos < (len(file_list) - 1):
        if filter_quality(file_list[pos + 3], args.qual) == True and \
                filter_bases(file_list[pos + 1], args.nbases) == True and \
                filter_short(file_list[pos + 1], args.len) == True:
            outputfile.write(file_list[pos] + '\n' + file_list[pos + 1] + '\n' + file_list[pos + 2] +
                             '\n' + file_list[pos + 3] + '\n')
            pos += 4
    outputfile.close()
    summary_file(args.sum_output, args.qual, args.nbases, args.len)


def main():
    parser = argparse.ArgumentParser(description="Reads Trimmer for FASTQ file")
    parser.add_argument("-in", help="fastq input file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="trimmed fastq filename", dest="output", type=str, required=True)
    parser.add_argument("-sum", help="summaryfilename", dest="sum_output", type=str,
                        default="Summaryfile")  # default or required true or not?
    parser.add_argument("-trim3", help="Number of nucleotides trimmed on 3' end", type=int, default=0)
    parser.add_argument("-trim5", help="Number of nucleotides trimmed on 5' end", type=int, default=0)
    parser.add_argument("-qual", help='Specifiy min. Quality of reads', type=int, default=0)
    parser.add_argument("-len", help="Minimum length of reads", type=int, default=0)
    parser.add_argument("-nbases", help="Maximum of unknown bases", type=int, default=1000)
    parser.set_defaults(func=run)
    args.func(args)


if __name__ == "__main__":
    main()
