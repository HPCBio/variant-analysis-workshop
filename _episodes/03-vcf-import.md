---
# Please do not edit this file directly; it is auto generated.
# Instead, please edit 03-vcf-import.md in _episodes_rmd/
source: Rmd
title: "Importing a VCF into Bioconductor"
teaching: 45
exercises: 15
questions:
- "How can I import and filter a VCF with Bioconductor?"
objectives:
- "Read genotypes and SNP metadata from a VCF into R."
- "Compress and index a VCF for quick access."
- "Use custom parameters to filter variants into a smaller file."
- "Choose which fields, samples, and genomic ranges to import into R."
keypoints:
- "Index the VCF file with indexTabix if you plan to only import certain ranges."
- "Use filterVcf to filter variants to a new file without importing data into R."
- "Use ScanVcfParam to specify which fields, samples, and genomic ranges you want to import."
---





~~~
library(GenomicRanges)
library(Rsamtools)
library(VariantAnnotation)
~~~
{: .language-r}

## Making the VCF file easy to access

We have a small VCF for this workshop, but they can be many gigabytes in size.
So, instead of reading straight from a plain text file, we will compress and
index the file.  This only needs to be done to the file once, even if you will
read data from the file many times in later R sessions.

To save time, the file you downloaded was zipped already.  If it wasn't, you
would zip it this way:


~~~
bg <- bgzip("data/hmp321_agpv4_chr1_subset.vcf")
~~~
{: .language-r}

> ## BGZIP format
>
> The BGZIP format is an extension of the GZIP format, so any command that
> would work on `.gz` files (like `zless`) works on `.bgz` files.  The
> BGZIP format was designed specifically for bioinformatics files to make
> them easy to index.  See http://www.htslib.org/doc/bgzip.html for more
> information.
{: .callout}

Instead we'll input the zipped file name directly.


~~~
bg <- "data/hmp321_agpv4_chr1_subset.vcf.bgz"
~~~
{: .language-r}

Now we'll build an index that will make it easy for R to read just particular
regions of the genome.

~~~
indexTabix(bg, format = "vcf")
~~~
{: .language-r}

~~~
[1] "data/hmp321_agpv4_chr1_subset.vcf.bgz.tbi"
~~~
{: .output}

Now we can make an object that points to the file for quick access, kind of like
the `FaFile` object we made in Episode 1.


~~~
myvcf <- VcfFile(bg)
~~~
{: .language-r}


{% include links.md %}