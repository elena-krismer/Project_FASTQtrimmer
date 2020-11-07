import sys
import statistics
import gzip
import re
import argparse

# user input - global or not  - bad for runtime
# nt3_input, nt5_input, threshold_reads_input, quality_input, n_bases_input = None, None, None, None, None

# counting triming

# filtered reads is this possible as global variable????
from pip._vendor.certifi.__main__ import args

filter_quality_count, filter_short_count, filter_bases_count = 0, 0, 0


def trim_user(seq_line, trim3, trim5):
    trim_line = seq_line[trim5:trim3]
    return trim_line


# quality score below 20 is considered low quality
# should user change
# the sliding window uses a relatively standard approach. this works by scanning from 5' end of the read and removes 3' end of the read when the a
# average quality of a group of bases drops below a specified threshold. This prevents a single weak base casusing the removal of subsequent high
# quality data,while still ensuring that ... (Bolger et al., 2014)
def trim_quality(seq_line, qual_line, phred):
    ascii_list = [ord(ascii_value) for ascii_value in qual_line]
    ascii_str = None
    if phred == '33':
        for val in ascii_list:
            val = int(val)
            while val in ascii_list > 53:
                ascii_str += ord(val)

    elif phred == '64':
        for val in ascii_list:
            val = int(val)
            while val in ascii_list > 84:
                ascii_str += ord(val)
    # cut the sequence the same length as the quality line
    return seq_line[0:len(ascii_str)], ascii_str


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
    n_number = (n_number + 1 for base == 'N' in base in seq_line)
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
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    args = parser.parse_args()

    typefile = re.search(r'\S*\.gz', args.input)
    if typefile:
        with gzip.open(args.file) as infile:
            file_list = [file_list.append(line.strip('\n')) for line in infile]
        infile.close()
    else:
        with open(args.file) as infile:
            file_list = [file_list.append(line.strip('\n')) for line in infile]
    infile.close()

    if file_list[0][0] != '@' or file_list[2][0] != '+' or re.search(r'[^ATGCN]', file_list[1]) is not None:
        print('Error in fileformat.')
        sys.exit(1)

    # trimming sequence
    pos_seq, pos_qual = 1, 3
    phred = detect_quality(file_list[6])
    while pos_qual < len(file_list):
        file_list[pos_seq] = trim_user(file_list[pos_seq], args.trim3, args.trim5)
        file_list[pos_qual] = trim_user(file_list[pos_qual], args.trim3, args.trim5)
        file_list[pos_seq] = trim_quality(file_list[pos_seq], file_list[pos_qual], phred)[
            0]  # first return from function
        file_list[pos_qual] = trim_quality(file_list[pos_seq], file_list[pos_qual], phred)[1]  # second return
        pos_seq += 4
        pos_qual += 4

    pos, trimmed = 0, 0
    outputfile = open(args.output)
    while pos < len(file_list):
        if filter_quality(file_list[pos + 3], args.qual) == True and \
                filter_bases(file_list[pos + 1], args.nbases) == True and \
                filter_short(file_list[pos + 1], args.len) == True:
            outputfile.write(file_list[pos] + '\n' + file_list[pos + 1] + '\n' + file_list[pos + 2] +
                             '\n' + file_list[pos + 3] + '\n')
        else:
            trimmed += 1
        pos += 4
    outputfile.close()
    summary_file(args.sum_output, args.qual, args.nbases, args.len)


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


def main():
    global func
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
    try:
        func = args.func
    except AttributeError:
        parser.error("too few arguments")
    func(args)
    # args.func(args)


if __name__ == "__main__":
    main()

