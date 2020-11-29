# Project FASTQ trimmer

The following program wil trim and filter your FASTQ file according to quality, length and unknown bases. The trimming based on quality, will trimm the ends of the read lower than a quality of 20. To run the programm you must a provide a FASTQ file in the standard FASTQ format (see Chapter [1](https://github.com/elena-krismer/Project_FASTQtrimmer/blob/master/Report/Report.md#1-Introduction)).

## Program Manual
This program allows you to filter and trim your FASTQ file. Additionally, a feature will provide you an overview of your FASTQ file, like average quality and number of unknown bases.



### 1. Trimming and Filtering 
The following program will trim and filter your FASTQ file according to quality, length and unknown (N's) bases. The trimming based on quality will trim the ends of the read lower than a quality of 20. To run the program you must provide a FASTQ file in the standard FASTQ format (see Chapter 5.3). Compressed as well as uncompressed files can be feed to the program. The output consists of two output files - a FASTQ file with filtered and trimmed reads and a summary file which contains information about the number of filtered + trimmed reads.

*Attention* :heavy_exclamation_mark: :warning:
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

##### For Filtering:

- **-qual** the minimum average quality of the read (default: Quality 20)
- **-length** the minimum length of the read 
- **-nbases** the minimum number of unknown bases

### 2. Statistics on FASTQ File
In case you are uncertain about setting the different parameters a statistics feature is implemented. This option will provide you with a statistic-summary file with information about the quality of the reads (average, the average quality of the worst and best 10% of the reads), the number of reads, the average length of the reads, and the total amount of the individual bases. Thus, with this information provided it will be easier to adjust parameters for trimming and filtering.

To perform statistics you must specify the name of your FASTQ input file (-in) and statistics-output file (-stat) and set the main-output file to false (-out False). Despite the commands '-in', '-out' and '-stat' not further commands should be used (see 5.3. Examples).

*Note*: You cannot perform trimming/filtering and statistics in one run.


### 3. Examples:

Following command trims 6 bases from each end of the read, filters all reads with a quality lower than 30, shorter than 50 nucleotides and more than two unknown bases.


```{p}
./fastqtrimmer.py -in S1.fastq -out S1_trim.fastq -sum SummaryS1_trim.txt -end3 6 -end5 6 -qual 30 -length 50 -nbases 2
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

