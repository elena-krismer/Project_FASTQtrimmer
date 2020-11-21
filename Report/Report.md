# Project
## 1. Introduction

The goal of this project is to generate a program, which trims Next-Generation Sequencing data.
## 2. Theory
## 3. Algorithm Design

This program relies on the uniform structure of FASTQ files: The first position the header followed by the sequence, the third line and the quality line. Thus, the file get read into a list and all following
operations were by calling these certain positions of the list (list position 1 for the sequence line, list position 3 for the quality line and so on).
In the first step, while iterating over the list with a while function, the trimming is performed and the quality scale determined.
After these modifications the reads are filtered. Whereby, the lines are passed to the functions which return Boolean Values. Only when
all three filter-functions (Mean Quality of the read, Number of unknown bases and the minimum length of the read) pass the test, thus return a True value, 
all four lines of the read get read into the output file. In case the read does not pass the test it will be not read into the outputfile
and counted as 'filtered'.

```{p}
main()
    read file into list
    while read list
        determine phred: position quality
        trim position 1(+4)/sequence and 3(+4)/quality
        trim quality 1(+4)/sequence and 3(+4)/quality
        count quality trims `

    while reading modified list
        filter quality, unknown bases, length if True:
            read first, second, third and fourth line into file
         else: count as filtered read
    write summary file with count of filtered and trimmed reads
```


## 4. Program Design
## 5. Program Manual`

## 6. Runtime Analysis

The main reasons for a slowdown in our runtime are the multiple function calls and
s. An alternative approach could be to store the lines in a Numpy Array. 
An alternative approach to improve runtime performance could be to use packages such as NumPy. Hereby, we would suggest storing
the lines in a NumPy Array and passing the NumPy Array(NumPy Array with Sequence string and NumPy Array with Quality string) to the functions, instead of passing each line separately to the function.

## 7. Discussion
## 8. References``````
