def filter(qual_line, quality, phred, threshold_reads):
    # converting to ascii values
    ascii_list = [ord(ascii_value) for ascii_value in qual_line]
    # sum quality and phred to get ascii value
    if statistics.mean(ascii_list) > (quality + phred) and len(qual_line) > threshold_reads:
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
        n_number = file_list[pos + 1].count('N')
        if filter(file_list[pos + 3], args.qual) == True and n_number < args.nbases:
            outputfile.write(file_list[pos] + '\n' + file_list[pos + 1] + '\n' + file_list[pos + 2] +
                             '\n' + file_list[pos + 3] + '\n')

        else:
            # counting filtered reads
            filtered += 1
        pos += 4
    outputfile.close()
