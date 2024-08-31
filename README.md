# Short Tandem Repeats (STRs) Procedure Notes

Last month I was given a copy of a friend's whole genome sequence (WGS) from Nebula Genomics, and tasked with looking for short tandem repeats (STRs) indicative/causative of diseases such as SBMA and Huntington's.

There is a lot of documentation, but much of it incomplete, some outdated, and rarely describing how to set up the data from the sequencing company. It took a month of long nights and dead ends, but now the process takes about a half-hour to set up, 11 hours of overnight crunching, and a half-hour to run and report.

And here it is for you, for what it's worth.

## Interpreting the results

The output will be a series of images, describing a known locus and the repeat counts, e.g. "AR n <= 14" indicates the androgen receptor locus, and that 14 repeats were found.

Refer to the STRipy database at https://stripy.org/database to interpret these results. Just go through each locus's report and compare the "n" to the database, see if any of them are in outside of the normal range and have a pathological length.
