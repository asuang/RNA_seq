#using hisat2, samtools and stringtie 
#hisat2 比对

ls -d */ |sed 's/\///' | uniq > filelist.txt
#ls *.fastq.gz  | sed 's/.sra.fastq.gz//' | uniq > filelist.txt
#ls *.fastq.gz  > filelist.txt

mkdir align
mkdir assembly
mkdir ballgown2
mkdir fastqc
species="C_elegans_Ensl_WBcel235"
cat filelist.txt | while read line
do 
	fastqc -o ./fastqc -t 16 ${line}/${line}_1.fq.gz
	fastqc -o ./fastqc -t 16 ${line}/${line}_2.fq.gz
done



cat filelist.txt | while read line
do
	echo $line	
	echo $line >> align/mapping_repo.txt
	#PE	
	hisat2 -p 16 --dta --rna-strandness RF -x /media/hp/disk1/song/Genomes/${species}/Sequence/hisat2/genome -1 ${line}/${line}_1.fq.gz -2 ${line}/${line}_2.fq.gz -S align/${line}.sam 2>> align/mapping_report.txt
	#SE	
	#hisat2 -p 16 --dta --rna-strandness RF -x /media/hp/disk1/song/Genomes/hg38/Sequence/hisat2/genome_tran -U ${line}.sra.fastq.gz -S align/$line.sam 2>> align/mapping_repo.txt
	#samtools 
	cd align
	samtools view -@ 16 -Sb ${line}.sam > ${line}.bam
	samtools sort -@ 16 ${line}.bam -o ${line}.sort.bam
	samtools index ${line}.sort.bam
	cd ../ballgown2	
	echo "### use stringtie to identify transcripts ###"
	mkdir $line
	stringtie -e -B -p 16 -G /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o $line/$line.FPKM.gtf ../align/${line}.sort.bam
	cd ..	
done


