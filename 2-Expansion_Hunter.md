# Setting up and Running Expansion Hunter

Having created our gene sequence file in BAM format, and deduplicated and indexed it, let's check for short tandem repeat (STR) diseases using Expansion Hunter.

Prerequisities:
- Follow the steps in document 1 to create your BAM file `ABC123PCW.sorted.dupmarked.bam`
- You will also be using the reference genome FASTA file `Homo_sapiens_assembly38.fasta` which you downloaded in that previous document.

This process was written for Linux systems, specifically Ubuntu 20.04 under Windows Subsystem for Linux (WSL).

Overview:
- Install Expansion Hunter, which exmaines genome files for STRs according to a list of known areas and repeats
- Install REViewer, which generates visualizations of the repeats found by EH
- Install GraphAlignmentViewer, which also generates visualizations of the repeats found by EH
- Run Expansion Hunter
- Run REViewer and GraphAlignmentViewer


## Installation

### Basic OS Stuff

```
sudo apt update
sudo apt upgrade

sudo apt install samtools
```

### Install Expansion Hunter

Website = https://github.com/Illumina/ExpansionHunter/

EH provides ready-to-run binaries which worked fine on WSL's Ubuntu 20.04 We just copy the binary into place where we can run it, and copy the variant catalog into our `/tmp/genome` so we can use it later.

```
# download and unpack
cd /tmp/genome
wget https://github.com/Illumina/ExpansionHunter/releases/download/v5.0.0/ExpansionHunter-v5.0.0-linux_x86_64.tar.gz
tar xvf ExpansionHunter-v5.0.0-linux_x86_64.tar.gz

# copy the binary and catalog files into place
sudo mv ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter ./ExpansionHunter
sudo chmod 755 ./ExpansionHunter
cp ExpansionHunter-v5.0.0-linux_x86_64/variant_catalog/hg38/variant_catalog.json ./ExpansionHunter_variantcatalog_hg38.json

# run it from CLI to see if it basically runs and shows usage info
./ExpansionHunter --help

# clean up the temp files
rm -fr ExpansionHunter-v5.0.0-linux_x86_64 ExpansionHunter-v5.0.0-linux_x86_64.tar.gz
```

### Install REViewer

Website = https://github.com/Illumina/REViewer

REViewer is prov ided as ready-to-run binaries, which worked fine on my WSL Ubuntu-20.04

```
# download and unpack
cd /tmp/genome
wget https://github.com/Illumina/REViewer/releases/download/v0.2.7/REViewer-v0.2.7-linux_x86_64.gz
gunzip REViewer-v0.2.7-linux_x86_64.gz

# copy it into place
mv REViewer-v0.2.7-linux_x86_64 ./REViewer
chmod 755 ./REViewer

# see if it runs and shows usage info
REViewer --help
```

### Install GraphAlignmentViewer

Website = https://github.com/Illumina/GraphAlignmentViewer

```
# GAV is written in Python and needs these other packages first
sudo apt install python3-matplotlib python3-pysam

# download via git
cd /tmp/genome
git clone https://github.com/Illumina/GraphAlignmentViewer.git

# make sure it runs and shows usage info
python3 ./GraphAlignmentViewer/GraphAlignmentViewer.py
```

## Setup and Run

```
# re-create the folder where results will go
cd /tmp/genome
rm -rf results
mkdir results

# run ExpansionHunter, generating a VCF report file and a BAM file of the examined areas
./ExpansionHunter \
    --reads ./ABC123PCW.sorted.dupmarked.bam \
    --reference ./Homo_sapiens_assembly38.fasta \
    --variant-catalog ./ExpansionHunter_variantcatalog_hg38.json \
    --output-prefix ./results/ehoutput

# sort and index the BAM file that EH generated, for the viewer tools 
samtools sort ./results/ehoutput_realigned.bam -o ./results/ehoutput_realigned_sorted.bam
samtools index ./results/ehoutput_realigned_sorted.bam

# run GraphAlignmentViewer, geneating a bunch of Results_GAV files for each examined locus
python3 ./GraphAlignmentViewer/GraphAlignmentViewer.py \
    --read_align ./results/ehoutput_realigned.bam \
    --gt_file ./results/ehoutput.vcf \
    --reference_fasta ./Homo_sapiens_assembly38.fasta \
    --variant_catalog ./ExpansionHunter_variantcatalog_hg38.json \
    --output_dir ./results/ --output_prefix Results_GAV

# run REViewer for each of the loci, saving a Results_REViewer file for each one
for locus in AFF2 AR ATN1 ATXN1 ATXN10 ATXN2 ATXN3 ATXN7 ATXN8OS C9ORF72 CACNA1A CBL CNBP CSTB DIP2B DMPK FMR1 FXN GIPC1 GLS HTT JPH3 NIPA1 NOP56 NOTCH2NL PABPN1 PHOX2B PPP2R2B RFC1 TBP TCF4; do
  echo $locus
  ./REViewer \
    --reads ./results/ehoutput_realigned_sorted.bam \
    --vcf ./results/ehoutput.vcf \
    --reference ./Homo_sapiens_assembly38.fasta \
    --catalog ./ExpansionHunter_variantcatalog_hg38.json \
    --locus $locus \
    --output-prefix results/Results_REViewer_{$locus}
done
```

You now have a `results/` folder full of images and reports, showing the STR lengths of everything that EH looked for.

## Diagnosis

Check out the various images (SVGs from REViewer and PNGs from GAV) against the database at https://stripy.org/database/ and see if any are in pathological ranges noted for each locus.
