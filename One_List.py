import sys
import statistics
import gzip
import re
import argparse
from pip._vendor.certifi.__main__ import args


def trim_user(seq_line, trim3, trim5):
    trim_line = seq_line[trim5:(trim3 + 1)]
    return trim_line

# trim quality according to quality of single base from 5'
def trim_quality(seq_line, qual_line, phred):
    qual_arr, qual_trim, count, pos, pos3 = bytearray(), bytearray(), 0, -1, 0
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
                pos3 += 1
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
                pos3 += 1
                pos -= 1

    # cut the sequence the same length as the quality line
    qual_dec = qual_trim.decode('utf-8')
    if pos3 != 0:
        qual_trim = qual_dec[0:-pos3]
        seq_trim = seq_line[:pos3]
        count = 1
    else:
        qual_trim = qual_dec

    if (len(seq_line) - len(qual_dec)) != 0:
        seq_trim = seq_line[(len(seq_line) - len(qual_dec)):]
        count = 1
    else:
        seq_trim = seq_line
    return seq_trim, qual_trim, count


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
    ascii_string = bytearray(b'ascii_string')
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


def run(args):
    file_list = list()
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    args = parser.parse_args()

    typefile = re.search(r'\S*\.gz', args.input)
    if typefile:
        with gzip.open(args.file) as infile:
            [file_list.append(line.strip('\n').replace(' ', '')) for line in infile]
        infile.close()
    else:
        with open(args.file) as infile:
            [file_list.append(line.strip('\n').replace(' ', '')) for line in infile]
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
                filter_bases_length(file_list[pos + 1], args.nbases, args.length) == True:
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
    parser.add_argument("-length", help="Minimum length of reads", type=int, default=0)
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

