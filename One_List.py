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



def trim_user(seq_line, trim3, trim5):
    trim_line = seq_line[trim5:trim3]
    return trim_line

# quality score below 20 is considered low quality
# the sliding window uses a relatively standard approach. this works by scanning from 5' end of the read and removes 3' end of the read when the a
# average quality of a group of bases drops below a specified threshold. This prevents a single weak base casusing the removal of subsequent high
# quality data,while still ensuring that ... (Bolger et al., 2014)


# trim quality according to quality of single base from 5'
def trim_quality(seq_line, qual_line, phred):
    ascii_list = [ord(ascii_chr) for ascii_chr in qual_line]
    ascii_str, count = None, 0
    if phred == '33':
        for ascii_val in ascii_list:
            ascii_val = int(ascii_val)
            while ascii_val in ascii_list > 53:
                ascii_str += chr(ascii_val)

    elif phred == '64':
        for ascii_val in ascii_list:
            ascii_val = int(ascii_val)
            while ascii_val in ascii_list > 84:
                ascii_str += chr(ascii_val)
    # cut the sequence the same length as the quality line
    seq_trim = seq_line[0:len(ascii_str)]
    # check whether sequence was trimmed for summary file
    if len(seq_line) != len(seq_trim):
        count = 1
    return seq_trim, ascii_str, count


# trim with moving window width of 4 bases
def moving_window(seq_line, qual_line, phred):
    # converting to values of character
    ascii_list = [ord(ascii_chr) for ascii_chr in qual_line]
    ascii_str, count, str, end = None, 0, 0, 3
    if phred == '33':
        while statistics.mean(ascii_list[str:end]) > 53:
                # chr() to convert to ascii_value back to character
                ascii_str.join([chr(ascii_value) for ascii_value in ascii_list[str:end]])
                str += 4
                end += 4
    else:
        while statistics.mean(ascii_list[str:end]) > 84:
                ascii_str.join([chr(ascii_value) for ascii_value in ascii_list[str:end]])
                str += 4
                end += 4
    # cut the sequence the same length as the quality line
    seq_trim = seq_line[0:len(ascii_str)]
    # check whether sequence was trimmed for summary file
    if len(seq_line) != len(seq_trim):
        count = 1
    return seq_trim, ascii_str, count


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


def filter_bases(seq_line, n_bases):
    n_number = 0
    n_number = (n_number + 1 for base == 'N' in base in seq_line)
    if n_number < n_bases:
        return True


def filter_short(seq_line, threshold_reads):
    if len(seq_line) > threshold_reads:
        return True


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
    pos_seq, pos_qual, trimmed = 1, 3, 0
    phred = detect_quality(file_list[6])
    while pos_qual < len(file_list):
        file_list[pos_seq] = trim_user(file_list[pos_seq], args.trim3, args.trim5)
        file_list[pos_qual] = trim_user(file_list[pos_qual], args.trim3, args.trim5)
        # saving trim_quality() return in variable
        func_return = trim_quality(file_list[pos_seq], file_list[pos_qual], phred)
        file_list[pos_seq], file_list[pos_qual] = func_return[0], func_return[1]
        trimmed += func_return[2]
        pos_seq += 4
        pos_qual += 4

    pos, filtered = 0, 0
    outputfile = open(args.output)
    while pos < len(file_list):
        if filter_quality(file_list[pos + 3], args.qual) == True and \
                filter_bases(file_list[pos + 1], args.nbases) == True and \
                filter_short(file_list[pos + 1], args.len) == True:
            outputfile.write(file_list[pos] + '\n' + file_list[pos + 1] + '\n' + file_list[pos + 2] +
                             '\n' + file_list[pos + 3] + '\n')
        else:
            # counting filtered reads
            filtered += 1
        pos += 4
    outputfile.close()

    with open('summaryfile', 'w') as sum_file:
        sum_file.write(
            "Total number of reads filtered: {0}. Reads trimmed {1}".format(filtered, trimmed))

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
    #args.func(args)


if __name__ == "__main__":
    main()

