import sys
import statistics
import gzip
import re
import argparse


def trim_user(seq_line, trim3, trim5):
    trim_line = seq_line[trim5:(trim3 + 1)]
    return trim_line


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


# functions trims quality lower than 20
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


def filter_quality(qual_arr, quality, phred):
    # sum quality and phred to get ascii value
    phred_ascii = quality + phred
    if len(qual_arr) != 0:
        if statistics.mean(qual_arr) > phred_ascii:
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

    qual_arr = bytearray()
    qual_arr.extend(map(ord, file_list[3]))
    phred = detect_quality(qual_arr)

    pos_seq, pos_qual, trimmed = 1, 3, 0
    phred = detect_quality(file_list[6])
    while pos_qual < len(file_list):
        file_list[pos_seq] = trim_user(file_list[pos_seq], args.trim3, args.trim5)
        file_list[pos_qual] = trim_user(file_list[pos_qual], args.trim3, args.trim5)
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

    pos, filtered = 0, 0
    with open('output.txt', 'w') as outputfile:
        while pos < (len(file_list) - 1):
            qual_arr = bytearray()
            qual_arr.extend(map(ord, file_list[pos + 3]))
            if filter_quality(qual_arr, args.qual, phred) == True and \
                    filter_bases_length(file_list[pos + 1], args.nbases, args.length) == True:
                outputfile.write(file_list[pos] + '\n' + file_list[pos + 1] + '\n' + file_list[pos + 2] +
                                 '\n' + file_list[pos + 3] + '\n')
            else:
                # counting filtered reads
                filtered += 1
            pos += 4
    outputfile.close()
