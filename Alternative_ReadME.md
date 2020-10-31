
# Project NGS - Code with lists
- Outputfile: gzip or 'normal' - call it **XXXX.trim.fastq.gz**

## Things defined by the user: - should these be global variables?
- triming on 3' and 5' end
- mean quality score 
- length of reads
- number of N bases 
### main()

    input: infile
reads all lines into a list, and performs the functions + reads the trimmed and filtered list into new file. By indexing the lines get passed to the function. 
So Heading[0], Seq-line[1], Third_line[2], Qual_line[3]


### def detect_quality():
The postion of the list containing the Q-score gets passed to the function and string gets converted into a list with corresponeding
ASCII values. The max and min value of the list is used to determine the range and so the Phred score.
    
    input: positon of the list (x+3) containing the quality 
    returns: phred scale (33/64)

### def trim_user():
    input: position of the list (x+3)(x+1) with sequence/quality
    output: trimmed string

trims the string from both sides according to users input


These trims from 5' can be done independently of each other in above order.
### def trim5_quality():
    input:
    output:

Filter out reads with a mean quality lower than specified after trimming.
### def filter_quality():
    input: positon of the list (x+3) containing the quality 
    output: *True* when the if-statements applies, else count has to be added to filter_quality_count (global variable)

Filter out reads that are shorter than specified after trimming.
### def filter_short():
    input: positon of the list *seq_line* *threshold_reads* user input
    output: *True* when if condition applies, else +1 triming variable 

Filter out reads that have a more than a specified number of N bases (unknown bases).
### def filter_bases():
    input:list 
    ouput: True or adding to global variable
 
 reads the N's in the given string 

until position 25 no errors -is our file filtered? in @ line Y=filtered,N otherwise
