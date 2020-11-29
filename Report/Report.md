# 22110 - Python and Unix for Bioinformaticians
<img align="right" width="300" height="100" src="DTU.jpg"><br /><br /><br /><br /> <br /><br /><br /><br /> <br /><br /><br />
# Project 8: Read trimmer for Next-Generation-Sequencing data<br /><br /><br /><br /> <br /><br /><br /><br /> <br /><br /><br />
# Students:
## Elena Krismer (s202425)
## Carolina Rocha (s203014)<br /> <br /><br /><br /><br /> <br /><br /><br /><br /> <br /><br /><br /><br /> <br /><br /><br />



# Content

### 1. [Introduction](#1)

### 2. [Theory](#2)

### 3. [Algorithm Design](#3)

### 4. [Program Design](#4)

### 5. [Program Manual](#5)

#### 5.1 [Main](#5.1)

#### 5.2 [Statistics](#5.2)
	
#### 5.3 [Examples](#5.3)


### 6. [Runtime Analysis](#6)

#### 6.1 [Big O](#6.1)

#### 6.2 [Visuale Profiling](#6.2)
	

### 7. [Conclusion](#7)

### 8. [Reference](#8)



To do/think of/ not forget?
- which time to we use when writting the theory future, present, conjunktiv; will be filter, should,... i am confused lol 'has to be adjusted to the ohred scle/will be adjust/must be adjusted/is adjusted???? ***After asking my metegnomics project friends and my mexican profrssor( no native reference lol) they said if is like a process is bettter present. I changed some lines there**
- yesterday was a merge conflict - we have to write the program design in pseudocode - i have added the pseudocode from yesterday, peter mentioned last time we shouldnt but the raw code in the report ***should we remove the gray code?** I think we can keep it because is very useful 

-yesterday I run the program with the peter file and I got a problem with the quality. I used another file I tested with the two programs and it runs (using another file) so I added into the directory, I''l move to the correct directory test

- should put the whole picture or a cutout in 6.2.? i am not sure if the whole picture is clear  *I think the whole picture explain the process*

- work on theory part


## 1. Introduction<a name="1">

Next-Generation Sequencing (NGS) has played and is still playing an important role in understanding biological mechanisms from a genomics perspective. In the early 2000s, the price to sequence a genome was immense. But with the decreasing sequencing costs, the genomic data production could be increased. Even though generating data became easier, computational storage and data analysis remain still a challenge. As the raw output genomic data contains sequencing errors, prep-processing is required to perform analysis. Different pipelines can be used to preprocess the data some of them share steps like a quality check, duplicate removal, and trimming reads. Read trimming is the process to remove low-quality bases or adapters while preserving the longest high-quality part of an NGS read. This step ameliorates mapping more reads to annotated genes, mitigates the effects of adapter contamination, plus it is widely assumed that trimming increases the accuracy of SNP calling and potentially could reduce the computational time (Didion et al., 2017; Del Fabbro et al., 2013; Bush, 2020). On the other hand, is the relevance of trimming in RNA-seq data is still discussed (Liao et Shi, 2020).

There have been several trimming tools developed. Given that, there is not one tool that simultaneously provides the accuracy, computational efficiency, and feature set to work with the types and volumes of data, the development and enhancement of trimming tools is still an emerging field. (Didion et al., 2017). The most common tools for trimming are Atropos, fastp, Trim Galore, and Trimmomatic (Bush, 2020).

There are two types of trimming based on 1) sequence and 2) quality. The first one can cut sequence adapters while the second one is based on the quality based on a Phred score. Both perspectives use a FASTQ file, which keeps the information of the sequencing and is conformed by: 

1. Header with the sequence identifier and information about the run and the cluster
2. The raw sequence (A,C,T, G and N)
3. "+" character separator sometimes followed by the header
4. Base quality score Phred +33 or +64 enconded, represented by ASCII characters


![](fastq.png)

*Figure 1-Structural example of a FASTQ format*


The fourth line in the read contains the quality score. The quality score (*Q*) decimals are logarithmically related to the error probability (*P*)( the probability that the base call is wrong) and is calculated as follow (Ochoa et al., 2013):

*Q = −10log10P*

The error probability, for each nucleotide, ranging from 0 to 1. Thus, 1 represents a probability of 100% for the nucleotide to be wrong and nucleotides with a p_error close to 0 to be correct( see *Table 1*). Bases with a high error probability are seen as an undetermined base and represented as 'N'.

The quality score is encrypted using the ASCII code into two systems, Phred +33 and +64. '33' and  '64' represent the first value in the scales, a quality score of 0 encoded as bytes (33 ASCII character = !; 64 ASCII character = @). The conversion between these two scales is relatively easy, as the quality score is encoded as decimals on Phred +64 scale, which is always 33 higher than the quality score encoded in decimals on the Phred +33 scale. For example, using the Phred +33 a quality of 20 will be represented by *“5”* which is the 53 number in ASCII code while *“T”* in +64 system (see the *Table 1*) (Ochoa et al., 2013).



![](qscores.gif)

*Table 1 Phred+33/+64 scale* - *source:usearchv11 page*
 

		
The quality, length, and the number of reads have a tremendous effect on the final results of experiments. Since the desired 'quality/quantity ratio' of the reads is depending on the further approach, we generated the program 'fastqtrimmer.py'.
		
This program allows to trim and filter Next-Generation Sequencing data from Illumina platforms. Whereby, the trimming and filtering parameters (quality, number of unknown bases and read length) can be defined by the user.




## Contribution

Report
Code


The final outputfile should only contain reads with a defined maximum of unknown bases, minimum average quality and minimum lenght of the sequence. The sequence/quality line must therefore meet all three criteria. When the filters are passed the four positions of the read are called and written into the ouputfile, else the read is counted as 'filtered' for the summaryfile.

## 2. Theory <a name="2">


As described in the introduction every read in a FASTQ file consists of four lines. This convention is the base of the program. Thus, the file gets read into a list and all following operations are performed by calling these certain positions of the list (list position 1 for the sequence line, list position 3 for the quality line).

The output of the program consists of two files the trimmed and filtered FASTQ file and the summary file, containing the count of trimmed and filtered reads. For trimming it has to be noticed that the position 'x' in the sequence line corresponds to position 'x' in the quality line. Thus when trimming the same amount of characters has to be trimmed from both lines. 

The final output file should only contain reads with a defined maximum of unknown bases, minimum average quality, and the minimum length of the sequence. The sequence/quality line must therefore meet all three criteria. When the filters are passed the four positions of the read will are called and written into the output file, else the read is counted as 'filtered' for the summary file.

Besides filtering and trimming the quality score has to be adjusted to the determined Phred scale. 

To not overwhelm the user with too many options, the trimming and quality parameters are optional.


 
 
## 3. Algorithm Design <a name="3">

For this program, a linear algorithm is used and the following programming structures are included:
- Sequences
- Binary Selection
- Repetition

The general idea of the algorithm is to transform the input file into a list, and while iterating in steps of four over the list performing several operations on the elements on the list. By using binary selection, the list elements are either written into the 'main' output file or counted for the secondary output file. The main algorithm and an example for a FASTQ read is represented in Figure *2*.

![](flowchart_update.png)

*Figure 2-Algorithm Scheme* 

*Note*: The binary selection between the trim/filter and the statistic operation is not mentioned.

## 4.Program Design <a name="4">

This program is written in Python.

### 4.1. <essentialMain part <a name="4.1">

##### Main steps:
To create a command-line interface the argparse library is used. To run the program the user must define the FASTQ filename and the name of the output file. All further commands are optional and the minimum quality is specified as 20.
After the arguments of the user got passed to the run-function, the following steps are conducted:

- **Reading into a list**: after reading the file the lines are stored into a list.

- **Determining Phred Scale**: input is the quality line as bytearray from 100th read.

```{p}
    detect_quality()
     	if mean(ASCII decimals) < 75:
      		phred_scale = 33
    	elif max(ASCII decimals ) >= 75:
      		phred_scale = 64
```

- **Trimming**: The 'main' trimming function (trimming_list) passes the strings to the function trim_user and trim_quality. The function trim_user slices the given number of characters from the input strings. Trim_quality takes the quality line (converted to a bytearray) and the sequence line as input.

```{p}
    def trimming_list()
    	while read list in an interval of four
         	convert ASCII characters to bytearray
         	trim_user: trimming bases list[position_sequence, position_quality]
         	trim_quality: list[position_sequence, position_quality]
         	count quality trims `
     
     def trim_user()
         slice characters from sequence and and bytearray
     
     def trim_quality()
         use determined phred scale to convert minimum quality score
         iterate over bytearray
            count number of characters to trim from 5' end
            count number of characters to trim from 3' end
         slice sequence and bytearray
         count trimming 
	 
    def trimming_list() return list to def_run()
    
```

- **Filtering and Writing in Outputfile**:

```{p}
    def write_outputfile()
   	 while reading trimmed list in interval of four:
          
	  def filter_quality
            	adjust quality score to phred scale, calculate mean of bytearray string
            	when average bigger than quality score return True
          
	  def filter_unkown_bases_length
            	count 'N', determine length
            	return True
         
         if all filters return True:
           	convert bytearray to ASCII string
           	write all four lines of read into outputfile
           
         else: count as filtered read
         
     def write_summaryfile()     
          write summary file with count of filtered and trimmed reads
```


**Following procedures should be mentioned explicitly, as they are fundamental for a valid output:**

- The determined Phred scale gets passed to the trim_quality and filter_quality function to adjust the quality score.

- The sequence and the quality line always get passed together to the trimming functions, to avoid a shift in the quality score of the bases.

- During the first iteration over the list, the quality line will be converted to a bytearray and will only be translated into an ASCII character string when the read gets written into the output file.



### 4.2. Statistics <a name="4.2">
	
Additionally to the trimming function, the program has a statistic function implemented. This operation provides instead of a trimmed and filtered FASTQ file a statistics-summary file for the given FASTQ file. Containing: the mean quality of a read, the mean quality of the best tenth and worst tenth of the reads, the average spot length, as well as the number of bases and the total number of reads.
This operation will only be conducted when it is explicitly specified by the user (see 5. Program Manual).

- **Statistics**

```{p}
    def run(takes argparse arguments)
    if outputfile false
       statistics
    
    def fastq_statistics()
      open file read into
        	sequence_list
        	quality_list
    
    def detect_quality() detect phred scale
    
    def statistics_numbases(sequence_list)
       		count bases in sequence_list
      		pass results to fastq_statistics
    
    def statistic_quality(quality_list - as bytearray)
         		create list with length of read
         		create list with average quality
       		sort average_quality_list
      		slice lower and upper tenpercent of the list and calculate mean
       		calculate mean with length_of_read_list
       		pass results to fastq_statistics
     
    def fastq_statistics()
    		write results from statistics_numbases and statistic_quality into summary	
       	
       
```



## 5. Program Manual <a name="5">

### 5.1. Trimming and Filtering <a name="5.1">
The following program will trim and filter your FASTQ file according to quality, length and unknown (N's) bases. To run the program you have to provide a FASTQ file in the standard FASTQ format (see Chapter 5.3). Compressed as well as uncompressed files can be fed to the program. The output consists of two output files - a FASTQ file with filtered and trimmed reads and a summary file which contains information about the number of filtered + trimmed reads.

*Attention* :heavy_exclamation_mark: :warning:
To make the script executable you must run following line:

```{p}
chmod +x fastqtrimmer.py
```

#### To run the program you have to specify:

- **-in** the input file 
- **-out** the name of the output file 

#### Further arguments are optional:

- **-sum** the name of the summary file (will give you a summary of the filtered and trimmed reads) (Defaultname: Summaryfile.txt)
- **-stat** ( *dont use this command for trimming* ) specifiying the name of your statistics output file (details in 5.2. Statistics on FASTQ File)

##### For Trimming:

- **-end5** the number of bases which should be trimmed from the 5´end
- **-end3** the number of bases which should be trimmed from the 3' end 

##### For Filtering:

- **-qual** the minimum average quality of the read (default: Quality 20)
- **-length** the minimum length of the read 
- **-nbases** the minimum number of unknown bases

### 5.2. Statistics on FASTQ File <a name="5.2">
In case you are uncertain about setting the different parameters, a statistics feature is implemented. This option will provide a statistic-summary file with information about the quality of the reads (average, the average quality of the worst and best 10% of the reads), the number of reads, the average length of the reads, and the total amount of the individual bases. Thus, with this information provided, it will be easier to adjust parameters for trimming and filtering.

To perform statistics you have to specify the name of your FASTQ input file (-in) and statistics-output file (-stat) and set the main-output file to false (-out False). Despite the commands '-in', '-out' and '-stat', no further commands should be used (see 5.3. Examples).

*Note*: You cannot perform trimming/filtering and statistics in one run.


### 5.3. Examples: <a name="5.3">

Following command trims 6 bases from each end of the read, filters all reads with a quality lower than 30, shorter than 50 nucleotides and more than two unknown bases.


```{p}
./fastqtrimmer.py -in S1.fastq -out S1_trim.fastq -sum SummaryS1_trim.txt -end3 6 -end5 6 -qual 30 -length 50 -nbases 2
```



In case you only want to trim the reads with a quality lower than 20 from each end and filter reads with a quality lower than 20, this command is enough:

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


## 6. Runtime Analysis <a name="6">

### 6.1. Big O <a name="6.1">
To evaluate the runtime in Big O terms a small overview of the functions, their complexity is defined in the following table. Most of the structures implemented shown linear time, which means that every single element from the input is visited exactly once, *O(n)* times. As the size of the input, N, grows our algorithm's run time scales exactly with the size of the input. N hereby indicates the number of lines/reads and m the length of the read.

| Function            | Big O analysis  |  
| ------------------- |     :---:    |
| Reading into a list | *O(n)*   | 
| detect_quality()  |   *O(1)*  | 
| trimming_list() | *O(n)* |
| trim_quality() | *O(2n)* reduced to  *O(n)* |
| filter_quality() | *O(n)* |
| filter_bases_length() | *O(n)* |
| write_outputfile() | *O(n)* |
| write_summary()   | *O(1)* |

*Table 2* Big O analysis in the main functions of the program. Where n is the number of lines, and m the length of the read.
*Note*: This table only represents an overview of the O-complexity.

O(n + 1 + n + 2n + n + n + n + 1) = *O(n)*

The length of a read, has a minor impact on the running time, thus can be ignored. The major factor for a linear increase in runtime is the number of lines (n)


### 6.2. Visual Profiling <a name="6.2">

To visualize the function calls and to get a better understanding for the runtime performance, we used the library [Python Call Graph](https://pycallgraph.readthedocs.io/en/master/). A cutout of those results is visible in Figure *3*. 

The trim and filter functions, do not distinguish significantly in runtime. The most complex operation, containing if/else statements and iteration is implemented in the trim_quality function and has the poorest runtime performance out of those four functions. As this is also the most complex function according to Big-O terms *O(2n)* (but still *O(n)*), an enhancement in this function will be beneficial. The 'cheapest' functions detetct_quality and write_summary will also remain that fast with bigger FASTQ files, as they have a Big-O complexity of *O(1)*.



![](pycallgraph_25.11.png)

*Figure 3-Cutout of the scheme generated by PyCallGraph. Script ran with a FASTQ file with 1000 reads.*



## 7. Discussion <a name="7">
	
One of the biggest limitations is that the program works using files from Illumina, it is not able to read files from 454, Nano, SOLID, or PacBio platforms due to quality detection.

The main bottleneck of the program is the detection of the Phred scale. The quality detection is extremely sensitive around the decimal 75 (= K), which is a quality score of 42 on Phred +33 scale and a quality score of 11 on Phred +64 scale. In case the read has considerably low quality (lower than 11) on a Phred scale +64, the Phred scale will be determined incorrectly as Phred scale +33. Since the quality of the first reads is commonly the lowest, we chose the quality of the 100th read (which is in a common FASTQ file still an early position) for detection. In further steps, there could be an error handling implemented, which uses the next read in case the quality scale of the first read can not be determined. As an alternative, another algorithm for the Phred scale determination should be considered. However, using the 100th position implies that a very small FASTQ file can not be fed to the program.

Further, a Big-O complexity of *O(n)* is considered a not-ideal, but desirable complexity. Instead of a linear increase of the runtime, a logarithmic increase should be aspired. However, we doubt that there is an applicable logarithmic algorithm design for this approach. 


The program is written in the dynamically-typed language Python, which means no explicit declaration of a variable is required but also implements that each variable contains extra information about their datatype. Thus, each element in our list contains its own information like the reference count and the datatype. As an alternative approach, a library such as NumPy could be used. As all lines in a FASTQ file are strings, storing them in a fixed-type array (NumPy Array) could increase efficiency. Hereby, we would suggest storing the lines in a NumPy Array and passing the NumPy Array (NumPy Array with thr sequence string and a NumPy Array with quality string) to the functions, instead of passing each line separately to the function (VanderPlas, 2016). This would also decrease the numerous function calls (visualized in 6.2.).

The main strength of the program is the easy handling. For every person, who knows how to use a command line.

Additionally, the modularization of the program allows changes(for instance only trimming or filtering of the file) without messing up the program.

Overall, the program is functional and provides the desired output.


## 8. References <a name="8">

Bush, S. J. (2020). Read trimming has minimal effect on bacterial SNP calling accuracy. *bioRxiv.*

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890.

Del Fabbro, C., Scalabrin, S., Morgante, M., & Giorgi, F. M. (2013). An extensive evaluation of read trimming effects on Illumina NGS data analysis. *PloS one*, 8(12), e85024.

Didion, J. P., Martin, M., & Collins, F. S. (2017). Atropos: specific, sensitive, and speedy trimming of sequencing reads. *PeerJ*, 5, e3720.

Liao, Y., & Shi, W. (2020). Read trimming is not required for mapping and quantification of RNA-seq reads at the gene level. *NAR Genomics and Bioinformatics*, 2(3), lqaa068.

usearch page: https://drive5.com/usearch/manual/quality_score.html

Ochoa I, Asnani H, Bharadia D, Chowdhury M, Weissman T, Yona G. (2013). QualComp: a new lossy compressor for quality scores based on rate distortion theory. *BMC Bioinformatics*. 2013;14:187.

VanderPlas, J. (2016). Python Data Science Handbook. *O'Reilly Media*.

## 9. List of Figures

1. Figure: Structural example of a Fastq format
2. Figure: Algorithm Scheme
3. Figure: Cutout of the scheme generated by PyCallGraph. Script run with a FASTQ file with 1000 reads

## 10. List of Tables

1. Table: Phred Scale
2. Table: Big-O Analysis
