 
mkdir "align"
mkdir "assembly"
mkdir "bigwig"

#比对
mkdir align
for i in `seq 1 12`
do
	hisat2 -p 16 --dta --rna-strandness RF -x /media/hp/disk1/song/pre_RNA_Seq/hg38_tran/genome_tran -1 G${i}_1.fq.gz -2 G${i}_2.fq.gz -S align/G${i}.sam
done

#安装samtools
wget https://sourceforge.net/projects/samtools/files/samtools/1.10/samtools-1.10.tar.bz2
tar -xjvf samtools-1.10.tar.bz
cd samtools-1.10
./configure --prefix=/media/hp/disk1/song
make
make install 
vi ~./bashrc
#加入
export PATH=/media/hp/disk1/song/bin:$PATH

#samtools 操作
cd align
for i in `seq 1 12`
do 
	samtools view -@ 16 -Sb G${i}.sam > G${i}.bam
	samtools sort -@ 16 G${i}.bam -o G${i}.sort.bam
	samtools index G${i}.sort.bam	
	rm G${i}.sam G${i}.bam
done
#特别注意rm那一步（确定能samtools了再用！！！不然会执行删除的）

#stringtie
cd ..
mkdir assembly
cd align
for i in `seq 1 12`
do 
	stringtie S${i}.sort.bam -o ../assembly/S${i}.gtf -p 16 -G /media/hp/disk1/song/RNA_Seq/hg38/Genes/genes.gtf
done

cd assembly
ls *.gtf > mergelist.txt
stringtie --merge -p 16 -G /media/hp/disk1/song/RNA_Seq/hg38/Genes/genes.gtf -o stringtie_merged.gtf mergelist.txt

#### compare the assembled transcripts to known transcripts  #####
gffcompare -r /media/hp/disk1/song/RNA_Seq/hg38/Genes/genes.gtf -o merged stringtie_merged.gtf
#或者-o 前面还可以加参数-G？

cd ..
ls *.fq.gz | sed 's/_[12].fq.gz//' | uniq > filelist.txt
mkdir ballgown

# make the gtf files with estimated RNA level # 
cat filelist.txt | while read LINE; do

	input=`echo "$LINE"`
	cd ballgown
 
	mkdir $input

	stringtie -e -B -p 8 -G ../assembly/stringtie_mergerd.gtf -o $input/$input.FPKM.gtf ../align/$input.sort.bam

	cd ..
	
done

cat filelist.txt | while read LINE; do

	input=`echo "$LINE"`
	cd ballgown2
 
	mkdir $input

	stringtie -e -B -p 8 -G /media/hp/disk1/song/raw_Data/hg38/Genes/genes.gtf -o $input/$input.FPKM.gtf ../align/$input.sort.bam

	cd ..
	
done

cat filelist.txt | while read LINE; do input=`echo $LINE`; cd ballgown2; mkdir $input; stringtie -e -B -p 16 -G ~/raw_Data/hg38/Genes/genes.gtf -o $input/$input.FPKM.gtf ../align/$input.sort.bam; cd ..; done





