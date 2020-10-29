**<font color="grey"><font size=40>UPSTREAM ANALYSIS of RNA-Seq </font></font>**
<font size=5><font color="grey"><p align="right">2020.10.27</p></font></font>
- [<font color="steelblue">Pipe for RNA-seq(polyA enrichment,mRNA) </font>](#-font-color--steelblue--pipe-for-rna-seq-ploya-enrichment-mrna----font-)
  - [<font size=4>1   quality control (fastQC)</font>](#-font-size-4-1---quality-control--fastqc---font-)
  - [<font size=4>2   mapping reads(RNA-seq mappers)</font>](#-font-size-4-2---mapping-reads-rna-seq-mappers---font-)
  - [<font size=4>3   sorting alignment and converting(samtools)</font>](#-font-size-4-3---sorting-alignment-and-converting-samtools---font-)
  - [<font size=4>4   assembly and quantitate(stringtie)</font>](#-font-size-4-4---assembly-and-quantitate-stringtie---font-)
  - [[RNA_seq.sh](https://github.com/asuang/RNA_seq/blob/master/RNA_seq.sh)]
- [<font color="steelblue">Pipe for RNA-seq(remove rRNA,total RNA) </font>](#-font-color--steelblue--pipe-for-rna-seq-remove-rrna-total-rna----font-)
  - [<font size=4>1   quality control (fastQC)</font>](#-font-size-4-1---quality-control--fastqc---font--1)
  - [<font size=4>2   mapping reads(RNA-seq mappers)</font>](#-font-size-4-2---mapping-reads-rna-seq-mappers---font--1)
  - [<font size=4>3   sorting alignment and converting(samtools)</font>](#-font-size-4-3---sorting-alignment-and-converting-samtools---font--1)
  - [<font size=4>4   assembly and quantitate(stringtie)</font>](#-font-size-4-4---assembly-and-quantitate-stringtie---font--1)
  - [<font size=4>5   compare (gffcompare)</font>](#-font-size-4-5---compare--gffcompare---font-)
  - [<font size=4>5  assembly and quantitate(stringtie)</font>](#-font-size-4-5--assembly-and-quantitate-stringtie---font-)
  - [[RNA_seq_total_RNA.sh](https://github.com/asuang/RNA_seq/blob/master/RNA_seq_total_RNA.sh)]
<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>


##  <font size=4>1   quality control (fastQC)</font>
Get the basic and quality information of the library.
```
fastqc -o ./fastqc -t 16 ${line}/${line}_1.fq.gz
```
##  <font size=4>2   mapping reads(RNA-seq mappers)</font>
Using the RNA-seq mappers , such as hisat2 , bowtie , bowtie2 or another , mapping the reads against the genome reference and identifying their genomic positions.
```
hisat2 -p 16 --dta --rna-strandness RF -x /media/hp/disk1/song/Genomes/${species}/Sequence/hisat2/genome -1 ${line}/${line}_1.fq.gz -2 ${line}/${line}_2.fq.gz -S align/${line}.sam 2>> align/mapping_report.txt
```
```
bowtie -p 16 -S /media/hp/disk1/song/Genomes/${species}/${species}_ref/bowtie/longest ${line}_20.fastq align/${line}.sam 2>> align/mapping_report.txt
```
```
bowtie2 -S ${input}.sam -p 16 -5 5 -3 10 -x /media/hp/disk1/song/Genomes/NC10/Sequences/WholeGenomeFasta/bowtie2/NC10 -1 ${input}_1.fq.gz -2 ${input}_2.fq.gz 2>> align/mapping_report.txt
```
##  <font size=4>3   sorting alignment and converting(samtools)</font>
Sorting the alignment by the genomic positions or names .
```
samtools view -@ 16 -Sb ${line}.sam > ${line}.bam
```
Converting .bam to .sam saves the storage.
```stringtie -e -B -p 16 -G /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o $line/$line.FPKM.gtf ../align/${line}.sort.bam
samtools sort -@ 16 ${line}.bam -o ${line}.sort.bam
```
Creating the index for the .bam.
```
samtools index ${line}.sort.bam
```
##  <font size=4>4   assembly and quantitate(stringtie)</font>
Reconstrcting all the isoforms from each genes and estimating their relative abundance.
```
stringtie -e -B -p 16 -G /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o $line/$line.FPKM.gtf ../align/${line}.sort.bam
```
## [RNA_seq.sh](https://github.com/asuang/RNA_seq/blob/master/RNA_seq.sh)

# <font color="steelblue">Pipe for RNA-seq(remove rRNA,total RNA) </font>
##  <font size=4>1   quality control (fastQC)</font>
 same with the step 1 described above
##  <font size=4>2   mapping reads(RNA-seq mappers)</font>
same with the step 2 described above
##  <font size=4>3   sorting alignment and converting(samtools)</font>
same with the step 3 described above
##  <font size=4>4   assembly and quantitate(stringtie)</font>
Assembly the new transcriptions are not in the reference.
```
stringtie -p 16 -G /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o $line.FPKM.gtf ../align/${line}.sort.bam
```
Merging every sample of the  .gtf to stringtie_merged.gtf.
```
stringtie --merge -p 16 -G /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o stringtie_merged.gtf mergelist.txt
```
##  <font size=4>5   compare (gffcompare)</font>
comparing the new .gtf with the reference .gtf.
```
gffcompare -r /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o merged stringtie_merged.gtf
```
##  <font size=4>5  assembly and quantitate(stringtie)</font>
Using the stringtie_merged.gtf as new reference.
```
stringtie -e -B -p 16 -G ../assembly/stringtie_merged.gtf -o $line/$line.FPKM.gtf ../align/${line}.sort.bam	
```

## [RNA_seq_total_RNA.sh](https://github.com/asuang/RNA_seq/blob/master/RNA_seq_total_RNA.sh)
