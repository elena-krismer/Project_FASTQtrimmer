#!/usr/bin/env python3
import sys
import gzip
import re
import argparse
from datetime import datetime
import timeit


# this functions is only for statistics
# counting number of bases
def statistics_numbases(seq_list):
    n_number, a_number, c_number, g_number, t_number = 0, 0, 0, 0, 0
    for item in seq_list:
        n_number += item.count('N')
        a_number += item.count('A')
        c_number += item.count('C')
        g_number += item.count('G')
        t_number += item.count('T')
    return n_number, a_number, c_number, g_number, t_number


# this function is only for statistics
# takes list of bytearray as input, calculating statistics
def statistic_quality(qual_list):
    avg_list, len_list = list(), list()
    for item in qual_list:
        len_list.append(len(item))
        avg_list.append(sum(item) / len(item))
    # sorting list to get lowest and highest tenth of quality
    avg_list.sort()
    length = len(avg_list)
    tenpercent = int(length * 0.1)
    avg_qual = round(sum(avg_list) / length)
    best_ten = round(sum(avg_list[-(tenpercent - 1):]) / len(avg_list[-(tenpercent - 1):]))
    worst_ten = round(sum(avg_list[:(tenpercent - 1)]) / len(avg_list[:(tenpercent - 1)]))
    entrie_length = round(sum(len_list) / len(len_list))
    return avg_qual, best_ten, worst_ten, entrie_length


# this function is only for statistics
# writing summary statistics file
def statistic_summary(statisticfile, infile, phred, avg_qual, best_ten, worst_ten, entries,
                      entrie_length, n_number, a_number, c_number, g_number, t_number):
    with open(statisticfile, 'w') as sum_file:
        sum_file.write(
            "# FASTQ Trimmer \n# Statistics FASTQ File \n# Python 3 \n# {0}\n# File: {1}\
             \nNumber of Reads: {2} \nAverage length of read: {3}\n\n## Number of bases: "
            "\n\tAdenin:\t{4}\n\tGuanin:\t{5}\n\tThymin:\t{6}\n\tCytosin:\t{7}\n\tUnknown Bases:\t{8}\n\n"
            "## Quality: \nPhred Scale: {9}\n\tAverage Quality of read:\t{10}"
            "\n\tAverage Quality of the best 10% of the reads:\t{11}"
            "\n\tAverage Quality of the worst 10% of the reads:\t{12}".format(
                datetime.now(), infile, entries, entrie_length, a_number, g_number, t_number,
                c_number, n_number, phred, avg_qual, best_ten, worst_ten))
    sum_file.close()


# this function is only for statistics
# main function for statistics
def fastq_statistics(infile, statisticfile):
    qual_list, seq_list = list(), list()
    typefile = re.search(r'\S*\.gz', infile)
    if typefile:
        try:
            infile = open(infile, mode='rt')
        except IOError as error:
            print('Can not open gzip file', str(error))
            sys.exit(1)
    else:
        try:
            infile = open(infile, 'r')
        except IOError as error:
            print('Can not open file', str(error))
            sys.exit(1)
    # reading into two lists, qual_list is list of bytearrays
    with infile:
        for x, line in enumerate(infile):
            if (x + 3) % 4 == 0:
                seq_list.append(line.strip('\n').replace(' ', ''))
            if (x + 1) % 4 == 0:
                line = line.strip('\n').replace(' ', '')
                qual_arr = bytearray()
                qual_arr.extend(map(ord, line))
                qual_list.append(qual_arr)
    infile.close()
    number_entries = len(seq_list)
    phred = detect_quality(qual_list[1])
    # applying phred scale to get quality
    func_return = statistic_quality(qual_list)
    avg_qual, best_ten, worst_ten = (func_return[0] - phred), (func_return[1] - phred), (func_return[2] - phred)
    entrie_length = func_return[3]
    func_return = statistics_numbases(seq_list)
    n_number, a_number, c_number = func_return[0], func_return[1], func_return[2]
    g_number, t_number = func_return[3], func_return[4]
    statistic_summary(statisticfile, infile, phred, avg_qual, best_ten, worst_ten, number_entries,
                      entrie_length, n_number, a_number, c_number, g_number, t_number)


# detection of phred scale using bytearray
# takes bytearray as input
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
    trim_line = seq_line[trim5:]
    if trim3 != 0:
        trim_line = trim_line[:-(trim3)]
    return trim_line


# functions trims quality lower than 20 from each end
def trim_quality(seq_line, qual_arr, phred):
    end5, count, pos, end3 = 0, 0, 0, 0
    # converting quality score from 20 to given phred scale
    phred_ascii = phred + 40
    # starting 5' end, counting characters to trim
    for pos in qual_arr:
        if pos < phred_ascii:
            pos += 1
            end5 += 1
        else:
            break
    # starting 3' end
    pos = len(qual_arr)
    for pos in qual_arr:
        if pos < phred_ascii:
            pos -= 1
            end3 += 1
        else:
            break

    # if qual_arr[pos] < phred_ascii:
    #   while qual_arr[pos] < phred_ascii:
    #      end5 += 1
    #      pos += 1
    # while pos >= 0:
    #     pos -= 1
    #     if qual_arr[pos] < phred_ascii:
    #         while qual_arr[pos] < phred_ascii:
    #             end3 += 1
    #             pos -= 1

    # trim 5' and 3' end, count for summaryfile
    if end3 != 0:
        qual_arr = qual_arr[:-end3]
        seq_line = seq_line[:-end3]
        count = 1
    if end5 != 0:
        qual_arr = qual_arr[end5:]
        seq_line = seq_line[end5:]
        count = 1
    return seq_line, qual_arr, count


# filter according to user specified mean quality - default 20
# takes bytearray as input
def filter_quality(qual_arr, quality, phred):
    # sum quality and phred to get ascii value
    phred_ascii = quality + phred
    if len(qual_arr) != 0:
        if (sum(qual_arr) / len(qual_arr)) > phred_ascii:
            return True


# filter according to specified length and min. number of unknown bases
def filter_bases_length(seq_line, n_bases, threshold_reads):
    n_number = seq_line.count('N')
    if n_number < n_bases and len(seq_line) > threshold_reads:
        return True


# main function for trimming
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
            # qual_arr = bytearray()
            # qual_arr.extend(map(ord, file_list[pos + 3]))
            if filter_quality(file_list[pos + 3], qual, phred) == True and \
                    filter_bases_length(file_list[pos + 1], nbases, length) == True:
                outputfile.write(file_list[pos] + '\n' + file_list[pos + 1] + '\n' + file_list[pos + 2] +
                                 '\n' + file_list[pos + 3].decode('utf-8') + '\n')
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
            "# FASTQ trimmer \n# Python 3 \n# {0}\n# File: {1}\
             \nTotal number of reads filtered: {2} \nTotal number of reads trimmed: {3}".format(
                datetime.now(), infile, filtered, trimmed))
    sum_file.close()


def run(args):
    file_list = list()
    infile, outputfile, summaryfile, statisticfile = args.input, args.output, args.sum_output, args.stat_output
    qual, trim3, trim5, length, nbases = args.qual, args.end3, args.end5, args.length, args.nbases
    # feature statistics
    if outputfile == 'False' and statisticfile != 'False':
        fastq_statistics(infile, statisticfile)
    # main function for trimming starts here
    else:
        typefile = re.search(r'\S*\.gz', infile)
        if typefile:
            try:
                infile = gzip.open(infile, mode='rt')
            except IOError as error:
                print('Can not open gzip file', str(error))
                sys.exit(1)
        else:
            try:
                infile = open(infile, 'r')
            except IOError as error:
                print('Can not open file', str(error))
                sys.exit(1)
        # reading file into list
        with infile:
            [file_list.append(line.strip('\n').replace(' ', '')) for line in infile]
        infile.close()
        # testing file format
        if file_list[0][0] != '@' or file_list[2][0] != '+' or re.search(r'[^ATGCN]', file_list[1]) is not None:
            print('Error in fileformat.')
            sys.exit(1)

        # determining quality scale
        qual_arr = bytearray()
        qual_arr.extend(map(ord, file_list[402]))
        phred = detect_quality(qual_arr)
        # trimming list
        func_return = trimming_list(file_list, trim3, trim5, phred)
        file_list = func_return[0]
        trimmed = func_return[1]
        # filter list, writing in file
        filtered = write_outputfile(file_list, outputfile, qual, phred, nbases, length)
        write_summary(summaryfile, trimmed, filtered, args.input)


def main():
    parser = argparse.ArgumentParser(description="Reads Trimmer for FASTQ file; "
                                                 "Feature: Performing statistics on FASTQ file ")
    parser.add_argument("-in", help="Filename for the input FASTQ file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="Filename for the trimmed FASTQ file", dest="output", type=str, required=True)
    parser.add_argument("-sum", help="Summaryfilename", dest="sum_output", type=str, default="Summaryfile.txt")
    parser.add_argument("-stat", help="Statistics for Input FASTQ file, define Filename. "
                                      "You must specify the outputfile as False (-out False) to perform statistics",
                        dest="stat_output", type=str, default="False")
    parser.add_argument("-end3", help="Number of bases trimmed on 3' end", type=int, default=0)
    parser.add_argument("-end5", help="Number of bases trimmed on 5' end", type=int, default=0)
    parser.add_argument("-qual", help='Specify min. quality of reads', type=int, default=20)
    parser.add_argument("-length", help="Minimum length of reads", type=int, default=0)
    parser.add_argument("-nbases", help="Maximum of unknown bases", type=int, default=1000)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
