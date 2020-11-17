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


def trim_quality(seq_line, qual_line, phred):
    qual_arr, qual_trim, count, pos, end3 = bytearray(), bytearray(), 0, -1, 0
    qual_arr.extend(map(ord, qual_line))
    if phred == 33:
        while pos < (len(qual_arr) - 1):
            pos += 1
            if qual_arr[pos] > 53:
                while pos < len(qual_arr):
                    qual_trim.append(qual_arr[pos])
                    pos += 1
        pos = (len(qual_trim) - 1)
        if qual_trim[pos] < 53:
            while qual_trim[pos] < 53:
                end3 += 1
                pos -= 1

    elif phred == 64:
        while pos < (len(qual_arr)):
            pos += 1
            if qual_arr[pos] > 84:
                while pos < len(qual_arr):
                    qual_trim.append(qual_arr[pos])
                    pos += 1
        pos = (len(qual_trim) - 1)
        if qual_trim[pos] < 84:
            while qual_trim[pos] < 84:
                end3 += 1
                pos -= 1

    # cut the sequence the same length as the quality line
    end5 = qual_trim.decode('utf-8')
    if end3 != 0:
        qual_trim = end5[0:-end3]
        seq_trim = seq_line[:end3]
        count = 1
    else:
        qual_trim = end5

    if (len(seq_line) - len(end5)) != 0:
        seq_trim = seq_line[(len(seq_line) - len(end5)):]
        count = 1
    else:
        seq_trim = seq_line
    return seq_trim, qual_trim, count

    # def main()
    # testing file format
    if file_list[0][0] != '@' or file_list[2][0] != '+' or re.search(r'[^ATGCN]', file_list[1]) is not None:
        print('Error in fileformat.')
        sys.exit(1)
