import argparse

def main():
    parser=argparse.ArgumentParser(description='Convert fast')
    parser.add_argument('-in', help='fasta ', dest='input', type=str)
    parser.add_argument('-out', help=' fastqoutput filename', dest='output')
    parser.add_argument('-qual', help='Quality score')
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

 if __name__ == '__main__':
        main()


