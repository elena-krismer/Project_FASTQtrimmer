import sys
import statistics
import gzip
import re
import argparse


def trim_user(seq_line, trim3, trim5):
    trim_line = seq_line[trim5:(trim3 + 1)]
    return trim_line


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


# functions trims quality lower than 20
def trim_quality(seq_line, qual_line, phred):
    qual_arr, qual_trim, count, pos, end3 = bytearray(), bytearray(), 0, -1, 0
    qual_arr.extend(map(ord, qual_line))
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


def filter_quality(qual_line, quality, phred):
    # converting to ascii values
    ascii_list = [ord(ascii_value) for ascii_value in qual_line]
    # sum quality and phred to get ascii value
    if statistics.mean(ascii_list) > (quality + phred):
        return True


def filter_bases_length(seq_line, n_bases, threshold_reads):
    n_number = seq_line.count('N')
    if n_number < n_bases and len(seq_line) > threshold_reads:
        return True

    # def main()
    # testing file format
    if file_list[0][0] != '@' or file_list[2][0] != '+' or re.search(r'[^ATGCN]', file_list[1]) is not None:
        print('Error in fileformat.')
        sys.exit(1)
