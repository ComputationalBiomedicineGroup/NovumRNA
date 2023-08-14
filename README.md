# NovumRNA: A non-canonical neoantigen prediction pipeline

## Happy to see you!

You want to predict non-canonical neoantigens (ncnas) from cancer RNA-seq data? You came to the right place! Thank you for having an interest in our tool NovumRNA! To ensure a smooth first experience with our tool, make sure to read this README carefully!

This README is meant to provide an overview and help you to run your first analysis using NovumRNA. For further information on parameters and available workflows, please refer to the user manual:
https://docs.google.com/document/d/1daVnVVYiOqdg7k4tWqshV0k8CYCQY0EBeGUC1gW4W8g/edit?usp=sharing

And if you are interested in how individual modules exactly work, take a look at the supplementary material:
TBA

Happy predicting! 

## Basic information about NovumRNA

NovumRNA takes RNA-seq FASTQ files from Tumor and optionally RNA-seq from healthy tissue to predict non-canoical neoantigens (ncnas).

In brief, NovumRNA takes single- or paired-end tumor RNA-seq FASTQ files as input. The pipeline also allows to re-run with intermediate files generated in a previous run. Reads are aligned to the reference genome with HISAT2 (Kim et al. 2019), or STAR (Dobin et al. 2013), to perform reference-guided transcript assembly using StringTie (Pertea et al. 2015). Tumor-specific and differentially expressed transcript fragments are identified by filtering against an internal normal control filtering database in the form of a capture BED file. Users can provide their own samples to extend the capture BED file. Transcript fragments need to fulfill expression and coverage criteria to be further considered, which can be defined by the user in the ```novumRNA.config``` file. Surviving transcript fragments are then translated and chopped into peptides of a specified length. Peptides are then filtered against the reference proteome. In parallel, the patient-specific human leukocyte antigen (HLA) class I and class II type is identified using OptiType (Szolek et al. 2014) and HLA-HD (Kawaguchi et al. 2017). Finally, tumor-specific peptides are tested for their binding affinity towards the patient’s HLA class I and class II molecules using netMHCpan-4.1 (Reynisson et al. 2020) and netMHCIIpan (Karosiene et al. 2013), to prioritize putative ncnas with higher likelihood of being presented to T cells. By providing additional annotation files, identified ncnas can be further categorized; NovumRNA by default provides a file to identify ncnas derived from ERVs.

## Prerequisites

NovumRNA is designed to run on a Linux based system.

NovumRNA is implemented as nextflow dsl2 pipeline which runs based on singularity containers, so you need to have nextflow >= 23.04.2.5870 and singularity >=  3.8.7 already installed on your system. 

If you don’t have them installed yet, the commands below may be used to install them.
 
```
curl -s https://get.nextflow.io | bash
curl -s https://get.singularity.io | bash
```
Or, if you have conda available, consider creating a dedicated Conda environment and install it there:

```
conda create -n novumRNA
conda activate novumRNA
conda install -c bioconda nextflow
conda install -c conda-forge singularity
```

## Installation

All the software is already installed within the singularity containers and is ready for you to use!

The only two things you need to do are:
1) Clone the repository to your system, which contains all the scripts you need.

```
git clone git@github.com:ausserh/NovumRNA.git
```

2) Download the NovumRNA resource bundle, containing the singularity containers, various reference files, testing data and many things more. The archive is 9 GB in size, when uncompressed it will increase to 16 GB, so make sure to have enough space.

```
wget TBA
tar -xzvf NovumRNA_resource_bundle.tar.gz
```
Did I say everything is already installed in the containers? Well, almost. Due to license restrictions, NovumRNA will install parts of its binding prediction module (IEDB tools) during its first run. You will need to accept a license, and then you are ready to go, NovumRNA takes care of the installation.

Something similar goes for the HLA-class II allele caller HLA-HD. We are not allowed to distribute its code. Therefore, HLA-HD needs to be already installed on your system. HLA-HD requires bowtie2, however, bowtie2 is already installed in the container.
So when you have HLA-HD installed, specify the path in the novumRNA.config file:
```HLAHD_DIR``` = "/path/to/hlhd.x.x.x/”

You don’t have HLA-HD or don’t want to use it? No big deal, the pipeline will run fine also without it, but you won’t get the HLA class II prediction. However, if you know which HLA-class II alleles you are interested in, there is also the possibility to provide them as input (see Usage: The csv samplesheet), the HLA-HD prediction will simply be skipped and you get your class II ncnas as well!

## Usage: The novumRNA.config file

Specifying the input, the output directory, changing references, aligner, or cut-offs, all of this happens in the ```novumRNA.config``` file, so this is your place to be!

More info on that in the manual, but here are the most essential positions you need to specify to get you started:

Whenever you specify a path, make sure to use the realpath and no symlink!


```input_ref``` = “path/to/resource/bundle” 
Specify here the path to your downloaded and unpacked resource bundle.

```novumrna``` = “path/to/novumra_git_repo” 
Specify here the path to your cloned NovumRNA git repository.

```outdir``` = “path/to/output_directory/” 
Specify here the path to where you want to have the pipeline output. If the directory is not present, it will be automatically created.

```input_fastq``` = “/path/to/input_samplesheet.csv” 
Specify here the path to your own input_samplesheet.csv (see next section).

However, ```input_fastq``` the config file is already set to:
"${params.input_ref}/samplesheet_CRC_fastq_sub.csv"

This points to the testing data present in the resource bundle. As a first run, we recommend using the testing data, so maybe keep this unchanged for the first run and then set it to your own samplesheet for further runs. The testing data consists of downsampled, colorectal cancer organoid derived RNA-seq data.

One more thing, depending on which scheduler you use, you need to define your scheduler parameters. This is done in so-called ```profiles```. Profiles are defined in the ```novumRNA.config``` file at the end. One is called singularity, you should leave this unchanged and always specify it for your runs (see Your first NovumRNA run). The profile called cluster defines how nextflow will submit the jobs, this depends on what scheduler you use. Modify these parameters: executor = 'slurm/sge/other' and clusterOptions = {#!/bin/bash … } or create your own profile.

## Usage: The csv samplesheet

A samplesheet is a format to specify your input data, nowadays widely used with nextflow pipelines. Your input files are given in form of a table, values separated by commas, hence a csv file.

The headers of this table consist of:
ID,Read1,Read2,HLA_types,HLA_types_II

* ID: Sample name, this ID will be used as output name
* Read1: forward read (FASTQ or gzipped FASTQ)
* Read2: reverse read (FASTQ or gzipped FASTQ, if left empty, single-end mode is used)
* HLA_types: A file containing already known HLA class I alleles, separated by commas, in this format: HLA-A*01:01,HLA-C*05:25 (see “Valid_HLAI_alleles.txt” in the repo bin).
If left empty, OptiType is run to predict the HLA class I alleles.
* HLA_types_II: A file containing already known HLA class II alleles, separated by commas, in this format: DPA1*01:03-DPB1*69:01,DRB1*15:13 (see “Valid_HLAII_alleles.txt” in the repo bin). If left empty, HLA-HD (if installed, like mentioned before) is run to predict the HLA class II alleles.

Every samplesheet needs to have these headers with these exact names. 

Here is how the samplesheet containing the testing data looks like:

```
ID,Read1,Read2,HLA_types,HLA_types_II
AK11_CRC01,AK11_CRC01_R1_combined_clean_rmdup.fastq,AK11_CRC01_R2_combined_clean_rmdup.fastq,,HLA_class_II_default_alleles.txt
```
As you can see, the column HLA_types can just be left empty, OptiType will run. For your test run, please specify the full path to these samples, all present in the resource bundle.

## Your first NovumRNA run

After you’ve set everything up in the config file and specified the paths in your own samplesheet, or the test samplesheet, run the pipeline from the command line like the following:

```
nextflow run path/to/repo/novumRNA.nf -profile cluster,singularity -w /path/to/work/  -c /path/to/repo/novumRNA.config -entry analysis
```
* ```profile```: cluster: your system specific scheduler options. singularity: options for singularity.
* ```-w```: Your working directory, where intermediate results will be stored
* ```-c```: Path to your novumRNA.config file
*```-entry```: use ```analysis``` for the standard prediction, see manual for other options.

Add “-resume” to your command, if something fails, and you need to re-run, the pipeline will take off at the last completed module, it will save you a lot of time!

As mentioned before, NovumRNA will install the IEDB toolkit on your first run.
After the first run, please specifiy ```IEDB_check``` in the novumRNA.config file to ${params.input_ref}/iedb/iedb_install_ok.chck. Like this the IEDB toolkit won’t be installed once more in a next run.

The default configurations in the novumRNA.config file will produce output:
hisat2 as aligner with a hisat2 index provided in the resource bundle (add ```--aligner star``` to the command to change to STAR as aligner, index will be created).
Class I ncnas of length 9 aas (add your length to ```peptide_length```  = "9 15" in the config).
Class II ncnas of length 15 aas, based on provided HLA class II alleles in the samplesheet, HLA-HD is not executed.

## Output explanation

All results can be found in the ```outdir``` directory you specified in the config. You can find results from each module there, the final list of ncnas is present in the file ```yourID_final_out_combined.tsv file```. One in the output folder ```Metadata_MHCI``` for class I and one in ```Metadata_MHCII``` for class II.

Only the IEDB installation, your own built indices (optional) and peptide references (optional) will appear in the resource bundle directory, more on this in the manual.

Output columns:

Note:
Transcript refers to the transcript where the ncna derives from
Exon refers to the exon within said transcript where the ncna derives from
TPM = Expression measurement, transcripts per million

* ```Chr``` = Chromosome where the ncna is found
* ```Peptide_START``` = Genomic start coordinate of ncna
* ```Peptide_STOP``` = Genomic stop coordinate of ncna
* ```Peptide``` = Ncna sequence
* ```Transcript``` = StringTie transcript ID
* ```STRAND``` = Genomic strand +/-
* ```ID``` = ID given in the samplesheet
* ```Mismatch``` = Nucleotide and position of SNPs in the ncna
* ```Annotation``` = Ncna class, either INTERGENIC, INTRON or Differential
* ```isoform_count``` = Number of transcript isoforms
* ```TPM_iso_fraction``` = (all isoform TPM sum / transcript TPM)*100
* ```Cov_within_VAF``` = Mean transcript exon coverage / ncna exon coverage
* ```E_START``` = Exon start coordinate
* ```E_STOP``` = Exon end coordinate
* ```NT_Overlap``` = Exon nucleotides overlap with capture BED region.
* ```Overlap_perc``` = Percentage of exon nucleotides overlap with capture BED region.
* ```E_Coverage``` = Exon read coverage reported by StringTie
* ```TPM``` = Expression of transcript in transcripts per million
* ```Transcript_ref``` = Translation reference, either “Protein” or “None”
* ```BAM_reads``` = Used reads (majority with same sequence) fully covering the ncna region in the BAM file
* ```BAM_reads_all``` = All reads (including differences, like SNPs) fully covering the ncna region in the BAM file
* ```Rank_HLA*``` = % Rank of ncna, binding prediction
* ```REF``` = Ncna aa sequence based on reference genome
* ```REF_NT``` = Ncna nt sequence based on reference genome
* ```NEW_NT``` = Ncna nt sequence including patient SNPs (if present)
* ```Annotation_2``` = provided extra annotation as BED file (default ERV annotation)
* ```NT_Overlap_2``` = Exon nucleotides overlap with extra annotation BED region

## Summary

* Make sure nextflow and singularity are installed
* Clone the NovumRNA repository
* Download and unpack the resource bundle
* In the config file, specify ```input_ref```, ```novumrna```, ```outdir``` and ```input_fastq```
* Set the paths to the testing files in the test samplesheet
* Adapt profile ```cluster``` to your system's scheduler
* Run with: nextflow run path/to/repo/novumRNA.nf -profile cluster,singularity -w /path/to/work/  -c /path/to/repo/novumRNA.config -entry analysis
* Set  “IEDB_check” in the config to path to “iedb_install_ok.chck”

## We hope it worked!
Following this small guide, we hope you have now a basic understanding how things work with NovumRNA and you should hopefully been able to run the tool using the provided test data, or your own data. If you are interested in using all the features of NovumRNA and learn more how to adapt it to your needs, read the official manual:
https://docs.google.com/document/d/1daVnVVYiOqdg7k4tWqshV0k8CYCQY0EBeGUC1gW4W8g/edit?usp=sharing
