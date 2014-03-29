========================================================

### CNVrd2: A package for measuring gene copy number, identifying SNPs tagging copy number variants, and detecting copy number polymorphic genomic regions


#### Install and use

Download the file 

> CNVrd2_1.1.4.tar.gz

Install the package

> R CMD INSTALL CNVrd2_1.1.4.tar.gz

Please see the file [**CNVrd2.pdf**](https://github.com/hoangtn/CNVrd2/blob/master/CNVrd2.pdf)

Window users can use the link of the Bioconductor Project:

http://www.bioconductor.org/packages/devel/bioc/html/CNVrd2.html


#### Notes: using the 1000 Genomes data

Please read information below or see the file [**using1000Genome.md**](https://github.com/hoangtn/CNVrd2/blob/master/using1000Genome.md) 



#### USING THE 1000 GENOMES DATA

This note describes some simple steps for using the data from the 1000 Genomes Project http://www.1000genomes.org/.

##   Bam files
Download an index file

```{}
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/alignment_indices/20130502.low_coverage.alignment.index
```

Obtain a list of bam files (the first column)

```{}
cat 20130502.low_coverage.alignment.index |awk '{print $1}'|grep '\.mapped' > listbam.txt
```

There are 2535 bam files on the page ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/. Therefore, we can make a list of these bam files and their full links.

```{}
awk '{print "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"$0}' < listbam.txt > listbamAndFullLinks.txt
```

We can choose a population (or multiple populations) to find tagSNPs. For example, here we choose the Mexican Ancestry in Los Angeles (MXL) population and find tagSNPs for FCGR3B gene.

```{}
cat listbamAndFullLinks.txt|grep "MXL" > listMXL.txt 
```
The gene is at chr1:161592986-161601753 http://www.ncbi.nlm.nih.gov/gene/2215, so we will use **samtools** (<a href="">Li et al. 2009</a> ) http://samtools.sourceforge.net/ to download a 1Mb region around the gene: chr1:161100000-162100000.

```{}
while read line
do
tempName=$(echo $line|awk -F"/" '{print $NF}')
samtools view -hb $line 1:161100000-162100000 > $tempName
done < listMXL.txt 
```

After downloading, we can use **samtools** to keep only reads mapped:

```{}
for file in $(ls *bam)
do

##Index bam file
samtools index $file

##Keep only read mapped
sammtols view -F 4 $file -b > temp.bam

mv temp.bam $file
done
```

A good website to understand SAM/BAM flags is http://picard.sourceforge.net/explain-flags.html

## VCF files

We can use **samtools** to obtain a vcf file (<a href="">Danecek et al. 2011</a> ) http://vcftools.sourceforge.net/ for all MXL samples above ***or*** we can use the SNP data of the 1000 Genomes Project.

Here, we use the data from the 1000 Genomes Project.

All VCF files can be obtained from this page http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/

We choose a region, for example chr1:161400000-161700000 flanking the FCGR3 gene to identify tagSNPs.

Download and index the vcf file by using "tabix" and "bgzip" command lines in the **tabix** tool (<a href="">Li, 2011</a> ) http://samtools.sourceforge.net/tabix.shtml:

```{}
tabix -h http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz 1:161400000-161700000  > chr1.161400000.161700000.vcf
```

Index the file:

```{}
bgzip chr1.161400000.161700000.vcf -c > chr1.161400000.161700000.vcf.gz

tabix -p vcf chr1.161400000.161700000.vcf.gz
```

## References

- P. Danecek, A. Auton, G. Abecasis, C.A. Albers, E. Banks, M.A. DePristo, R.E. Handsaker, G. Lunter, G.T. Marth, S.T. Sherry,  others,   (2011) The variant call format and VCFtools.  <em>Bioinformatics</em>  <strong>27</strong>  (15)   2156-2158
- H. Li, B. Handsaker, A. Wysoker, T. Fennell, J. Ruan, N. Homer, G. Marth, G. Abecasis, R. Durbin,  others,   (2009) The sequence alignment/map format and SAMtools.  <em>Bioinformatics</em>  <strong>25</strong>  (16)   2078-2079
- Heng Li,   (2011) Tabix: fast retrieval of sequence features from generic TAB-delimited files.  <em>Bioinformatics</em>  <strong>27</strong>  (5)   718-719
