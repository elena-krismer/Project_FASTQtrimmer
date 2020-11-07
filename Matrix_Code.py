import sys
import statistics
import gzip
import re
import argparse
import numpy as np
import time

# user input - global or not  - bad for runtime
# nt3_input, nt5_input, threshold_reads_input, quality_input, n_bases_input = None, None, None, None, None

# counting triming

# filtered reads is this possible as global variable????
from pip._vendor.certifi.__main__ import args

filter_quality_count, filter_short_count, filter_bases_count = 0, 0, 0


def trim_user(seq_list, trim3, trim5):
    seq_list = [elem[trim5:trim3] for elem in seq_list]
    return seq_list


def main():
    file_list = list()
    trim3 = 50
    trim5 = 10
    with open('BRISCOE_0069_BD18RUACXX_L2_1_pf.fastq.txt', 'r') as infile:
        for line in infile:
            file_list.append(line)
        infile.close()
    # trimming sequence
    length = len(file_list)

    file_matrix = [file_list[i:i + 4] for i in range(0, len(file_list), 4)]
    t_matrix = [[file_matrix[j][i] for j in range(len(file_matrix))] for i in range(len(file_matrix[0]))]
    t_matrix[1] = trim_user(t_matrix[1], trim3, trim5)
    t_matrix[3] = trim_user(t_matrix[3], trim3, trim5)
    print(t_matrix[1])

if __name__ == "__main__":
    main()
    import timeit
    print(timeit.timeit("main()"))