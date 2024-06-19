# Preparing the Data from the Sequencer

## Overview

This first stage is to create the complete BAM file (the uncompressed, realigned gene sequence) from the FASTQ files (which are the sequenced data with repeats and quality control information). This BAM file is the ready-to-analyze sequence data, having been de-duplicated, aligned to a reference genome, and sorted and indexed.

Our input files will be:
- The pair of FASTQ files, gzip-compressed, provided from Nebula Genomics, e.g. ABC123PCW_1.fq.gz and ABC123PCW_2.fq.gz
- A reference genome file fior the HG38 genome, in FASTA format

Tools used will be:
- minimap2, a tool for this one step, reconstructing tghe sequence data (in SAM format) from the FASTQ files
- samtools, a tool for converting, indexing, and otherwise handling gene sequences in SAM/CRAM/BAM format

You will also need:
- 400 GB of available disk space

The operating system will be:
- Linux, in this case Ubuntu 20.04; Windows Subsystem for Linux (WSL) works fine too

## Data Prep

Working directory will be `/tmp/genome`

```
mkdir /tmp/genome
```

Move the two FASTQ files into position:
```
cd /tmp/genome
mv ~/Downloads/ABC123PCW_1.fq.gz .
mv ~/Downloads/ABC123PCW_2.fq.gz .
```

Download the reference genome FASTA and its index. See https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0 for the full listing.
```
cd /tmp/genome
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
```

### Install Tools

The tools we need are in the Ubuntu repositories, so this is easy:

```
sudo apt install -y samtools minimap2
```

## Convert

This uses minimap2 to merge the FASTA files according to the reference genome, into a gene file in SAM format. However, we don't need the interm ediate SAM file and samtools are great for piping data from one step to another, to improve processing speed and minimize the use of (very large) temporary files.

```
cd /tmp/genome
rm -f ABC123PCW*.bam ABC123PCW*.bai _temp*.bam

minimap2 -t 4 -a -x sr ./Homo_sapiens_assembly38.fasta ABC123PCW_1.fq.gz ABC123PCW_2.fq.gz | \
  samtools fixmate -O bam,nthreads=3 - - | \
  samtools sort -T _temp -O bam,nthreads=6 - | \
  samtools markdup -O bam,nthreads=6,level=9 - ABC123PCW.sorted.dupmarked.bam

samtools index ABC123PCW.sorted.dupmarked.bam
```

The above commands assume that you have 8 logical CPU cores. If you have more or less you can change the `nthreads=` and `-t 4` to give those steps more/fewer processing threads.

The conversion pipeline took 10 hours on my reasonably-priced eight-core i7 system, and the indexing took about 10 minutes. Your times will vary.


## Done

The resulting BAM file `ABC123PCW.sorted.dupmarked.bam` is now suitable for use with Expansion Hunter, STRipy, or others.
