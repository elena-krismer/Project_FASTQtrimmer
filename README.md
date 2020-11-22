# Project_NGS_1.0.0
Python_unix

# Link to google colab

https://drive.google.com/drive/folders/1PBD5KWx19O-2vwSVXtcV-2Kg2oZr60Da?usp=sharing 

Following program wil trim and filter your FASTQ file according to quality, length and unknown bases. 

## To run the program you must specify:

- **-in** the inputfile 
- **-out** the outputfile 

### Further arguments are optional but are adviced to be specified:

- **-sum** the name of the summaryfile (will give you a summary of the filtered and trimmed reads)
#### For trimming:

- **-end5** the number of bases which should be trimmed from the 5´end
- **-end3** the number of bases which should be trimmed from the 3' end 

#### For Filter:

- **-qual** The minimum average quality of the read
- **-length** the minimum length of the read 
- **-nbases** the minimum of unknown bases

## Addention!
To make the script executable you must run following line:

```{p}
chmod +x fastqtrimmer.py
```


