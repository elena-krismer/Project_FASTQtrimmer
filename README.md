# Project_NGS_1.0.0
Python_unix

# Link to google colab

https://drive.google.com/drive/folders/1PBD5KWx19O-2vwSVXtcV-2Kg2oZr60Da?usp=sharing 

Following program wil trimm and filter your FASTQ file according to quality, length and unknown bases. The trimming based on quality, will trimm the ends of the read lower than a quality of 20. To run the programm you must a provide a FASTQ file in the standard FASTQ format (see Chapter X)

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

## Addention!
To make the script executable you must run following line:

```{p}
chmod +x fastqtrimmer.py
```


