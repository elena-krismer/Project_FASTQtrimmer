import sys
import statistics
import gzip
import re

# user input - global or not  - bad for runtime
nt3_input, nt5_input, threshold_reads_input, quality_input, n_bases_input = None, None, None, None, None

# counting triming

# filtered reads is this possible as global variable????
filter_quality_count, filter_short_count, filter_bases_count = 0, 0, 0


def trim_user(seq_line):
    trim_line = seq_line[:-nt3_input]
    trim_line = trim_line[nt5_input:]
    return trim_line


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
    if n_number > n_bases:
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


def main():
    file_list = list()
    with open('BRISCOE_0069_BD18RUACXX_L2_1_pf.fastq.txt', 'r') as infile:
        for line in infile:
            file_list.append(line)
        infile.close()
    # trimming sequence
    pos_seq = 1
    while pos_seq < len(file_list):
        file_list[pos_seq] = trim_user(file_list[pos_seq])
        pos_seq += 4

    # trimming quality scale
    pos_qual = 3
    while pos_qual < (len(file_list) - 1):  # -1 or not?
        file_list[pos_qual] = trim_user(file_list[pos_qual])
        pos_qual += 4
    q_scale = detect_quality(file_list[6])

    pos = 0
    while pos < (len(file_list) - 1):
        if filter_quality(file_list[pos + 3], quality_input) == True and \
                filter_bases(file_list[pos + 1], n_bases_input) == True and \
                filter_short(file_list[pos + 1], threshold_reads_input) == True:
            finalfile.write(file_list[pos] + '\n' + file_list[pos + 1] + '\n' + file_list[pos + 2] + \
                            '\n' + file_list[pos + 3] + '\n')
            pos += 4


def summary_file(summaryfile):
    with open('summaryfile', 'w') as sum_file:
        global filter_quality_count, filter_bases_count, filter_short_count
        filtered = filter_quality_count + filter_bases_count + filter_short_count
        sum_file.write("Total number of reads filtered: " + filtered + '\n' + \
                       filter_quality_count + "reads with a low quality than" + quality_input + '\n' \
                       + filter_bases_count + 'reads with more than' + filter_bases_count + 'bases')


if __name__ == "__main__":
    main()
