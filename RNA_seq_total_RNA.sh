#using hisat2, samtools and stringtie 
#hisat2 比对
species="C_elegans_Ensl_WBcel235"

#ls -d */ |sed 's/\///' | uniq > filelist.txt
#ls *.fastq.gz  | sed 's/.sra.fastq.gz//' | uniq > filelist.txt
#ls *.fastq.gz  > filelist.txt
mkdir assembly
mkdir ballgown2_total

cat filelist.txt | while read line
do
	echo ${line}
	cd assembly
	stringtie -p 16 -G /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o $line.FPKM.gtf ../align/${line}.sort.bam
	cd ..
done

cd assembly
ls *.gtf > mergelist.txt
stringtie --merge -p 16 -G /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o stringtie_merged.gtf mergelist.txt
gffcompare -r /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o merged stringtie_merged.gtf

cd ../ballgown2_total
cat ../filelist.txt | while read line
do	
	echo "### use stringtie to identify transcripts ###"
	mkdir $line
	stringtie -e -B -p 16 -G ../assembly/stringtie_merged.gtf -o $line/$line.FPKM.gtf ../align/${line}.sort.bam	
done


