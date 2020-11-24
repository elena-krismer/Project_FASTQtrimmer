# Project
## 1. Introduction

The goal of this project is to generate a program, which trims Next-Generation Sequencing data based on quality. 

- purpose of our programm 

## 2. Theory
Next Generation Sequencing has played an important role to understand the biology mechanisms under a genomics perspective.In the earlies X the price of sequence a genome was  very high but with time, the sequencing cost has decreased and the genomic data production has increased.Generating data became easier but not the computational storage and data analysis. This output genomic data is raw and contains error sequencing in order to perform analysis downstream it must be pre-processed. There are different pipelines that could be used to preprocess the data some of them share steps like quality check, duplicated removal,  and  trimming reads.Read trimming is the process to remove low quality bases or adapters while preserving the longest high quality part of a NGS read. Trimming step led to more reads mapping to annotated genes, mitigate the effects of adapter contamination, widely assumed to increase the accuracy of SNP calling and potentially could  reduce the computational time(Didion,J.P *et al*., 2017; Del Fabbro et al., 2013;  J. Bush S., 2020) on another hand there are studies where still discussing the trimming effect in RNA-seq data suggesting that read trimming is a redundant process in the quantification of RNA-seq expression data(Liao Y and Shi W., 2020). 

Didion and colleagues mention that several trimming tools had been developed  however there is not one that simultaneously provides the accuracy, computational efficiency and feature set  to work with the types and volumes of data (Didion,J.P *et al*., 2017) reason why different tools are still emerging. The most common tools for trimming are Atropos, fastp, Trim Galore, and Trimmomatic(Bush S., 2020).

There are two types of trimming based on 1) sequence and 2) quality. The first one is able to cut sequence adapters while the second one nucleotides based on the quality based on a prhed score. Both perspectives use a fastq file, this file keeps the information of the sequencing and is conformed by: 

1. Header with the sequence identifier and information about the run and the cluster
2. The sequence (A,C,T, G and N)
3. "+" character separator sometimes followed by the header
4. Base quality score Phred +33 or +64 enconded, represented by ASCII characters


![](fastq.png)

*Figure 1-Structural example of a Fastq format*


The quality score is encrypted using the ascii code into two systems phred 33 and 64. The first one adds the 33 into the quality, the second one works in the same way so instead of adding 33 you add 64. For example using the phred+33 a quality of 20 will be represented by *“5”* who is the 53 number in ASCII code while *“T”* in +64 system (see the *Table 1*)


![](qscores.gif)

*Table 1 Phred+33/+64 scale*
 
 Every ASCII character represents the error propability of each nucleotide to be correct, the values are from 0 to 1,  as lower the value more certain that the nucleotide is correct while 1 means that the base is certainly wrong ( see *Table 1*). This values closer to 1 in the nucleotide sequence appears like an undeterminated based (represented as 'N'). In order to remove the reads with a specified number of N's and trim low quality nucleotides at 3' and 5' we present this program.
 
 
 
## 3. Algorithm Design


![](BasicAlgorithm_BetterQuality.png)

*Figure 2 Algorithm Scheme*

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
ADD the diagram




has to be keeped in mind: different phred scales, structure of a fastq file, simulatenous trimming of quality and sequence line 

## 5. Program Manual

### 5.1. Trimming and Filtering 
Following program will trimm and filter your FASTQ file according to quality, length and unknown (N's) bases. The trimming based on quality, will trimm the ends of the read lower than a quality of 20. To run the programm you must a provide a FASTQ file in the standard FASTQ format (see Chapter X). The output consist of two outputfiles - a fastq file with filtered and trimmed reads and a summaryfile which contains information about the number of filtered + trimmed reads.


*Attention* :heavy_exclamation_mark:
To make the script executable you must run following line:

```{p}
chmod +x fastqtrimmer.py
```

#### To run the program you must specify:

- **-in** the inputfile 
- **-out** the name of the outputfile 

#### Further arguments are optional:

- **-sum** the name of the summaryfile (will give you a summary of the filtered and trimmed reads) (Defaultname: Summaryfile.txt)
- **-stat**  *dont use this command for trimming*  specifiying the name of your statistics output file (details in 5.2. Statistics on FASTQ File)

##### For Trimming:

- **-end5** the number of bases which should be trimmed from the 5´end
- **-end3** the number of bases which should be trimmed from the 3' end 

##### For Filter:

- **-qual** the minimum average quality of the read (default: Quality 20)
- **-length** the minimum length of the read 
- **-nbases** the minimum number of unknown bases

### 5.2. Statistics on FASTQ File
In case you are uncertain about setting the different parameteres a statistics feature is implemented. This option will provide you with a statisitc-summaryfile with information about the quality of the reads (average, average quality of the worst and best 10% of the reads), number of reads, the average length of the reads and the total amount of the individual bases. Thus, with this information provided it will be easier to adjust parameters for trimming and filtering.

To perform statistics you must specify the name of your fastq inputfile (-in) and statistics-outputfile (-stat) and set the main-outputfile to false (-out False). Despite the commands '-in', '-out' and '-stat' not further commands should be used (see 5.3. Examples).

*Note*: You cannot perform trimming/filtering and statistics in one run.


### 5.3. Examples:

Following command trims 6 bases from each end of the read, filters all reads with a quality lower than 30, shorter than 50 nucleotides and more than two unknown bases.


```{p}
./fastqtrimmer.py -in Sample1.fastq -out Sample1_trimmed.fastq -sum SummarySample1_trimming.txt -end3 6 -end5 6 -qual 30 -length 50 -nbases 2
```



In case you only want to trimm the reads with a quality lower than 20 from each end and filter reads with a quality lower than 20, this command is enough:

```{p}
./fastqtrimmer.py -in Sample1.fastq -out Sample1_trimmed.fastq
```


Performing statistics on a FASTQ file:
```{p}
./fastqtrimmer.py -in Sample1.fastq -out False -stat Sample1_Statistics.txt
```


To get an overview over the commands you can use, use following command:
```{p}
./fastqtrimmer.py -h
```


## 6. Runtime Analysis

The main reasons for a slowdown in our runtime are the multiple function calls and
s. An alternative approach could be to store the lines in a Numpy Array. 
An alternative approach to improve runtime performance could be to use packages such as NumPy. Hereby, we would suggest storing
the lines in a NumPy Array and passing the NumPy Array(NumPy Array with Sequence string and NumPy Array with Quality string) to the functions, instead of passing each line separately to the function.

To visualize the function calls and get a better understanding for the runtime performance we used the library [Python Call Graph](https://pycallgraph.readthedocs.io/en/master/). A cutout of those results are visible in Figure X. Since the programm passes the strings individually to the functions, the amount of function calls is noticably high. Considering that for creating this scheme a fastq-file with 1000 reads was used and a common fastq file is much bigger, an reduction of these function calls should be strived. Further there is a significant difference between the several trim and filter functions in runtime. Especially, the functions trim_quality() and filter_quality() require a huge amount of time to compute. Comparing those two functions, which both get a bytearray as input and perform stat

![](overview_runtime.png)

*Figure X* Cutout of the scheme generated by PyCallGraph. Script run with a FASTQ file with 1000 reads. 

Figure X. depict




## 7. Discussion


Runtime, a lot of function calls, finding alternative to list list of lists, tuples, arrays,... COMPARE WITH EXISTING ALGORITHMS

Relying on uniform format of fastq file 

Quality of single residue trimming is 20 cant be changed by user, in further approach make it optional
We were surprised by the huge amount of runtime the filter_quality function required


Quality of single residue trimming is 20 cant be changed by user, in further approach make it optional

## 8. References
Bush, S. J. (2020). Read trimming has minimal effect on bacterial SNP calling accuracy. *bioRxiv.*
Del Fabbro, C., Scalabrin, S., Morgante, M., & Giorgi, F. M. (2013). An extensive evaluation of read trimming effects on Illumina NGS data analysis. *PloS one*, 8(12), e85024.
Didion, J. P., Martin, M., & Collins, F. S. (2017). Atropos: specific, sensitive, and speedy trimming of sequencing reads. *PeerJ*, 5, e3720.
Liao, Y., & Shi, W. (2020). Read trimming is not required for mapping and quantification of RNA-seq reads at the gene level. *NAR Genomics and Bioinformatics*, 2(3), lqaa068.

## 9. List of Figures

1. Figure: 
2. Figure:
