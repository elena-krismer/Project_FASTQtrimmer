# Project_NGS_1.0.0

The following program wil trimm and filter your FASTQ file according to quality, length and unknown bases. The trimming based on quality, will trimm the ends of the read lower than a quality of 20. To run the programm you must a provide a FASTQ file in the standard FASTQ format (see Chapter [1](https://github.com/elena-krismer/Project_FASTQtrimmer/blob/master/Report/Report.md#Project)).

## To run the program you must specify:

- **-in** the inputfile 
- **-out** the name of the outputfile 

### Further arguments are optional but are adviced to be specified:

- **-sum** the name of the summaryfile (will give you a summary of the filtered and trimmed reads) (Defaultname: Summaryfile)

#### For Trimming:

- **-end5** the number of bases which should be trimmed from the 5´end
- **-end3** the number of bases which should be trimmed from the 3' end 

#### For Filter:

- **-qual** The minimum average quality of the read (default: Quality 20)
- **-length** the minimum length of the read 
- **-nbases** the minimum of unknown bases

## Attention!!
To make the script executable you must run following line:

```{p}
chmod +x fastqtrimmer.py
```

### Examples:

Following command trims 6 bases from each end of the read, filters all reads with a quality lower than 30, shorter than 50 nucleotides and more than two unknown bases.


```{p}
./fastqtrimmer.py -in Sample1.fastq -out Sample1_trimmed.fastq -sum SummarySample1_trimming.txt -end3 6 -end5 6 -qual 30 -length 50 -nbases 2
```


In case you only want to trimm the reads with a quality lower than 20 from each end, use following command:

```{p}
./fastqtrimmer.py -in Sample1.fastq -out Sample1_trimmed.fastq
```

