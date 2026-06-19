# SR_benchmarking
This routine is designed to assist in the analysis of short-read platform error rates by determining the total number of single-nucleotide, insertion, and deletion mismatches that occur after mapping with Donor-specific Assemblies (DSAs). Since DSAs are built to be idealized representations of a specific sample's reference, any mismatches found between the mapped reads and DSA itself can be thought to be a result of errors in the sequencing process. 

# Table of Contents
- [Pre-requisites](#pre-requisites)
- [Usage](#usage)
- [Outputs](#outputs)
- [Base quality filtering](#base-quality-filtering)
- [Low-complexity region analysis](#low-complexity-region-analysis)
- [Logfile](#logile)
  - [logfile example](#logfile-example)
- [Creating the appropriate BAM file](#creating-the-appropriate-BAM-file)

## Pre-requisites

### Programs
The SR_benchmarking script requires you to have `samtools` installed on your machine so that it can be run from within the Perl script with the command `samtools view ...`

**[samtools](https://github.com/samtools/samtools)** (v1.20 or higher) (https://github.com/samtools/samtools/blob/develop/INSTALL)

If you wish to perform an analysis on only the mismatches that fall into the low-complexity as defined by the [GA4GH database](https://github.com/usnistgov/giab-stratifications) consortium, you will also need to have `bedtools` installed.

**[bedtools](https://github.com/arq5x/bedtools2)** (v2.31.1 or higher)


### BAM Formatting
The SR_benchmarking script is designed to be run on a BAM file that contains both a CIGAR string (standard) and the optional [MD tag string](https://davetang.org/wiki/tiki-index.php?page=SAM) which contains details on the mismatching positions. 

If you do not already have a BAM file in this format, see **Creating the appropriate BAM file** from a set of fastq files below.


# Usage
```bash
perl study.mismatches_by_quality.pl -prefix [unique_run_identifier] -output_directory [path_to_output_files] -bam [BAM_file] -dump_details [optional; default=no]
```


# Outputs
Multiple output files are created in the `-output_directory` containing details of each mismatch type and a summary of the number of instances of each mismatch type. This is where the base-quality information is kept (Note: any "X" in the flanks column in the single_nucleotide_details file means the mismatch occurs at the first or last base of the read):

Single-nucleotide mismatch (SNM) files:
- **single-nucleotide.details.[run_prefix].bed** (location and details for each SNM)
- **single-nucleotide.summary.[run_prefix].txt** 

Insertion/deletion (indel) mismatch files:
- **insertion.details.[run_prefix].bed** (location and details for each insertion mismatch)
- **deletion.details.[run_prefix].bed** (location and details for each deletion mismatch)
- **indel.summary.[run_prefix].txt** 

  Reported base quality (BQ) for indels:
  - insertions: average BQ of all inserted bases
  - deletions: average BQ of the two bases flanking the deletion (use with caution, if at all)

Other summary files:
- **raw_mismatch_rates.summary.[run_prefix].txt** (summary of the raw mismatch rates -- no BQ filtering)
- **insert_size_distribution.[run_prefix].txt**
- **read_length_distribution.[run_prefix].txt**
- **mean_and_median.read_length_and_insert_size.[run_prefix].txt** (summary of RLs and ISs)
- **total_adjusted_read_length.[run_prefix].txt** (total bps analyzed -- used for the denominator in mismatch calculations)

# Base quality filtering
The SNM and indel mismatch details output bed files are designed to be easily filtered by the user to perform any desired **base-quality dependent mismatch rate** calculations.

For example, to determine the number of single-nucleotide mismatches that have a BQ>=30:
```bash
grep -v quality [path_to_output_files]/single-nucleotide.details.[run_prefix].bed | cut -f 4 | awk '{if ($0 >= 30) SUM++ } END{print SUM}'
```
The resulting value divided by the **total adjusted read length** is the **mismatch rate**.

# Low-complexity region analysis

A **low-complexity region analysis** can be performed by filtering the SNM and indel mismatch details output bed files on the appropriate genomic low-complexity regions found in the [GA4GH database](https://github.com/usnistgov/giab-stratifications). Low-complexity region bed files for the HG002 and COLO829BL Donor-specific Assemblies (DSAs), as well as for the GRCh38 reference will be provided in the Git Large File Storage (Git LFS) resources directory soon.

For example, to determine the number of insertion mismatches that fall into the low-complexity regions of the reference your sample was mapped to, first use `bedtools` to create a new file containing only events that fall into the defined low-complexity regions:
```bash
bedtools intersect -header -wa -a [path_to_output_files]/insertion.details.[run_prefix].bed] -b [corresponding_reference_low-complexity_bed] > [path_to_output_files]/low_complexity.insertion.details.[run_prefix].bed]
```

Then you can use the **base quality filtering** step from above to do any further analysis of those events based on their BQ values.

# Logfile
**logfile.[run_prefix].txt**

If `-dump_details yes` is set at the command line, the logfile will contain details for each read with an indel (sample below) -- which should always contain a few concurrent single-nucleotide mismatches as well.


However, the end of the logfile should always contain a summary of the mismatch rates.

### logfile example
>>>>>>>>>>>>>>>>>>>>>>>
locale:	chr10_MATERNAL:33-140  
readID:	LH00266:7:227CWMLT3:4:1122:41877:25017  
 CIGAR:	49M6I51M45S  
 MATCH:	36C8A4A11A9C5C8A12  
  MapQ:	6  
  
  read:	CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCATAACCCTACCCCTCCCAATCACCCTAACC  
 quals:	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF-FFFF-FFF---F---F-F-----55-5--5-555  
 index:	123456789|123456789|123456789|123456789|123456789|123456789|123456  
  
event_type&emsp;	pos&emsp;	ref&emsp;	alt&emsp;	flanks&emsp;	BQ&emsp;	(finalQ	leftOrRight)  
mismatch&emsp;	37&emsp;	C&emsp;	A&emsp;	C_T&emsp;	12  
mismatch&emsp;	46&emsp;	A&emsp;	C&emsp;	A_C&emsp;	12  
insertion&emsp;	50&emsp;	_&emsp;	TCCCAA&emsp;	C_C&emsp;	37,12,12,12,12,12_or_12,37,12,12,12,12&emsp;	16&emsp;	TAKING RIGHT-JUSTIFIED QSUM  
mismatch&emsp;	57&emsp;	A&emsp;	C&emsp;	T_A&emsp;	20  


# Creating the appropriate BAM file
This is best accomplished by using **[bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)** with the following command:
```bash
bwa-mem2 mem -t 4 [reference_fasta] [read_1_fastq] [[read_2_fastq]] > [output_SAM_file]
samtools sort -@ 4 -O BAM -o [output_BAM_file] [output_SAM_file]
samtools index -@ 4 [output_BAM_file]
```
