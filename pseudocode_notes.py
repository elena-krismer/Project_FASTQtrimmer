import argparse
import re
import gzip
parser = argparse.ArgumentParser()
parser.add_argument('file')
args = parser.parse_args()



def open_file(filename):
   typefile = re.search(r'\S*\.gz', args.file)
   if typefile:
      with gzip.open(args.file,'r') as f:
         lines=f.readlines()
         head=[item[:-1] for item in lines[::4]]
         read=[item[:-1] for item in lines[1::4]]
         qual=[item[:-1] for item in lines[3::4]]

         a=list(zip(head, read,qual))
         #alist=list(a)
         #print(a)

   else:
      with open(args.file) as f:
         lines=f.readlines()
         head=[item[:-1] for item in lines[::4]]
         read=[item[:-1] for item in lines[1::4]]
         qual=[item[:-1] for item in lines[3::4]]
         
         ai=list(zip(head, read,qual)) # what if we save the lines in a tuple? they will looks like this : 
#('@HISEQ_HU01:89:H7YRLADXX:1:1101:8895:6882 1:N:0:ATCACG', 'AAGGAGATGTGGGCGGGGAGAGGACGGGGGTCAGAAGACCGAGGGCGACCTCGAGGCGAGGGCGGGACAGCGGCGGGGGGTGAAGTACGCATGCGGATTCC', '
#@@CFFFFFHHGHHJIJJIEEEFFFFDDDDD:@BCDDCACDDDBBBDDDDDDDBDBBDDBDDDDDDDBBDDBDDBBD9>&500(:(+4(+&+(((5&0&)(+')



##### this part should be adapt but it seems like works
n_bases = 1  # maybe the user can specify it??
min_len = 10  # maybe the user can specify it??
seq_line = ["ATGNATGTGTGTGTGTGAGTAGTGAGTATATGTGTANNN", "NNNGTGTGTGTATGATAGTATGAGTAG", "GTGTGTTAGTAGTGTGAGTAG",
            "ATGGTAGTAGTACCCGTACGTA", "NNAGGGTAGTAGTA"]  ### Line added to probe with a list of reads
for member in seq_line:
    n_number = member.count('N')
    #print(n_number)  ## just to see if it works
    if n_number < n_bases:  # maybe here we can add defaoutl <= 0 just to keep the reads wihout N's??
        # print(item)
        if len(member) > min_len:
            print(member)

