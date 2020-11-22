#!/usr/bin/env python3
import sys
import statistics
import gzip
import re
import argparse
from datetime import datetime


# this is the code with additonal functions

# detection of phred scale using bytearray
def detect_quality(ascii_string):
    # phred 33 range= 33-75
    if max(ascii_string) <= 75 and min(ascii_string) < 59:
        phred_scale = 33
    # phred 64 range = 64 - 106
    elif max(ascii_string) > 75 and min(ascii_string) >= 64:
        phred_scale = 64
    else:
        print("Error in determining quality scale")
        sys.exit(1)
    return phred_scale


# trimming user specified 3' and 5' end
def trim_user(seq_line, trim3, trim5):
    trim_line = seq_line[trim5:-(trim3 + 1)]
    return trim_line


# functions trims quality lower than 20 from each end
def trim_quality(seq_line, qual_arr, phred):
    qual_trim, count, pos, end3 = bytearray(), 0, -1, 0
    # converting quality score from 20 to given phred scale
    phred_ascii = phred + 20
    # starting add 5' end reading quality bytearray to new bytearray
    while pos < (len(qual_arr) - 1):
        pos += 1
        if qual_arr[pos] > phred_ascii:
            while pos < (len(qual_arr) - 1):
                qual_trim.append(qual_arr[pos])
                pos += 1
    # starting add 3' end determining numbers of characters to trim
    pos = (len(qual_trim) - 1)
    while pos >= 0:
        pos -= 1
        if qual_trim[pos] < phred_ascii:
            while qual_trim[pos] < phred_ascii:
                end3 += 1
                pos -= 1
    # cut the sequence the same length as the quality line
    untrimmed_len = len(seq_line)
    qual_trim = qual_trim.decode('utf-8')
    # trim number of characters 3' end seq_line and qual_line
    if end3 != 0:
        qual_trim = qual_trim[0:-end3]
        seq_line = seq_line[:end3]
        count = 1
    # trim 5' end from seq_line
    if (untrimmed_len - len(qual_trim)) != 0:
        seq_line = seq_line[(untrimmed_len - len(qual_trim)):]
        count = 1
    return seq_line, qual_trim, count


# filter according to user specified mean quality - default 20
def filter_quality(qual_arr, quality, phred):
    # sum quality and phred to get ascii value
    phred_ascii = quality + phred
    if len(qual_arr) != 0:
        if statistics.mean(qual_arr) > phred_ascii:
            return True


# filter according to specified length and min. number of unknown bases
def filter_bases_length(seq_line, n_bases, threshold_reads):
    n_number = seq_line.count('N')
    if n_number < n_bases and len(seq_line) > threshold_reads:
        return True


def trimming_list(file_list, trim3, trim5, phred):
    pos_seq, pos_qual, trimmed = 1, 3, 0
    while pos_qual < len(file_list):
        file_list[pos_seq] = trim_user(file_list[pos_seq], trim3, trim5)
        file_list[pos_qual] = trim_user(file_list[pos_qual], trim3, trim5)
        # saving trim_quality() return in variable
        qual_arr = bytearray()
        qual_arr.extend(map(ord, file_list[pos_qual]))
        func_return = trim_quality(file_list[pos_seq], qual_arr, phred)
        # check whether the whole string got stringed, delete whole read in that case
        if func_return[0] is not None:
            file_list[pos_seq], file_list[pos_qual] = func_return[0], func_return[1]
        else:
            file_list[(pos_seq - 1):(pos_seq + 2)] = []
        trimmed += func_return[2]
        pos_seq += 4
        pos_qual += 4
    return file_list, trimmed


def write_outputfile(file_list, outputfile, qual, phred, nbases, length):
    pos, filtered = 0, 0
    with open(outputfile, 'w') as outputfile:
        while pos < (len(file_list) - 1):
            qual_arr = bytearray()
            qual_arr.extend(map(ord, file_list[pos + 3]))
            if filter_quality(qual_arr, qual, phred) == True and \
                    filter_bases_length(file_list[pos + 1], nbases, length) == True:
                outputfile.write(file_list[pos] + '\n' + file_list[pos + 1] + '\n' + file_list[pos + 2] +
                                 '\n' + file_list[pos + 3] + '\n')
            else:
                # counting filtered reads
                filtered += 1
            pos += 4
    outputfile.close()
    return filtered


# creating summaryfile with number of trimmed and filtered reads
def write_summary(summaryfile, trimmed, filtered, infile):
    with open(summaryfile, 'w') as sum_file:
        sum_file.write(
            "# Next Generation Sequencing Trimming \n# Python 3 \n# {0}\n# File: {1}\
             \nTotal number of reads filtered: {2} \nTotal number of reads trimmed: {3}".format(
                datetime.now(), infile, filtered, trimmed))
    sum_file.close()


def run(args):
    file_list = list()
    infile, outputfile, summaryfile = args.input, args.output, args.sum_output
    qual, trim3, trim5, length, nbases = args.qual, args.end3, args.end5, args.length, args.nbases

    # reading file into list
    typefile = re.search(r'\S*\.gz', infile)
    if typefile:
        try:
            with gzip.open(infile, mode='rt') as infile:
                [file_list.append(line.strip('\n').replace(' ', '')) for line in infile]
            infile.close()
        except IOError as error:
            print('Can not open gzip file', str(error))
            sys.exit(1)
    else:
        try:
            with open(infile) as infile:
                [file_list.append(line.strip('\n').replace(' ', '')) for line in infile]
            infile.close()
        except IOError as error:
            print('Can not open file', str(error))
            sys.exit(1)

    # testing file format
    if file_list[0][0] != '@' or file_list[2][0] != '+' or re.search(r'[^ATGCN]', file_list[1]) is not None:
        print('Error in fileformat.')
        sys.exit(1)

    # determining quality scale
    qual_arr = bytearray()
    qual_arr.extend(map(ord, file_list[6]))
    phred = detect_quality(qual_arr)

    # trimming list
    func_return = trimming_list(file_list, trim3, trim5, phred)
    file_list = func_return[0]
    trimmed = func_return[1]

    # filter list, writing in file
    filtered = write_outputfile(file_list, outputfile, qual, phred, nbases, length)
    write_summary(summaryfile, trimmed, filtered, args.input)


def main():
    parser = argparse.ArgumentParser(description="Reads Trimmer for FASTQ file")
    parser.add_argument("-in", help="fastq input file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="trimmed fastq filename", dest="output", type=str, required=True)
    parser.add_argument("-sum", help="Summaryfilename", dest="sum_output", type=str, default="Summaryfile.txt")
    parser.add_argument("-end3", help="Number of bases trimmed on 3' end", type=int, default=0)
    parser.add_argument("-end5", help="Number of bases trimmed on 5' end", type=int, default=0)
    parser.add_argument("-qual", help='Specifiy min. Quality of reads', type=int, default=20)
    parser.add_argument("-length", help="Minimum length of reads", type=int, default=0)
    parser.add_argument("-nbases", help="Maximum of unknown bases", type=int, default=1000)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
