## unzip files
cd $1
gunzip *.gz
mv CC*_1.fq $1.1.fq
mv CC*_2.fq $1.2.fq

## mapping PE-seq MNase-seq data and filtering reads
bwa mem -t 10 /mnt/Storage/data/BWA/mm9.fa ./$1.1.fq ./$1.2.fq > $1.sam
samtools view -bht /mnt/Storage2/home/wangchenfei/annotations/mm9/chromInfo_mm9.txt $1.sam > $1.bam
bamtools filter -mapQuality ">=10" -in $1.bam -out $1.Q10.bam

## Tag extension and pileup, 37bp
python bin/NEPTUNE_quality_control.py $1.Q10.bam mm9
python bin/NEPTUNE_quality_control.py $1.Q10.uniq.bed 37
sort -k1,1 -k2,2n $1.Q10.uniq.nucleo.37bp.bed > $1.Q10.uniq.nucleo.37bp.tmp
mv $1.Q10.uniq.nucleo.37bp.tmp $1.Q10.uniq.nucleo.37bp.bed
genomeCoverageBed -i $1.Q10.uniq.nucleo.37bp.bed -g /mnt/Storage2/home/wangchenfei/annotations/mm9/mm9.len -bga -scale $2 > $1.Q10.nucleo.37bp.bdg
bedGraphToBigWig $1.Q10.nucleo.37bp.bdg /mnt/Storage2/home/wangchenfei/annotations/mm9/mm9.len $1.Q10.nucleo.37bp.bw

## Tag extension and pileup, 73bp
python bin/NEPTUNE_quality_control.py $1.Q10.uniq.bed 73
sort -k1,1 -k2,2n $1.Q10.uniq.nucleo.73bp.bed > $1.Q10.uniq.nucleo.73bp.tmp
mv $1.Q10.uniq.nucleo.73bp.tmp $1.Q10.uniq.nucleo.73bp.bed
genomeCoverageBed -i $1.Q10.uniq.nucleo.73bp.bed -g /mnt/Storage2/home/wangchenfei/annotations/mm9/mm9.len -bga -scale $2 > $1.Q10.nucleo.73bp.bdg
bedGraphToBigWig $1.Q10.nucleo.73bp.bdg /mnt/Storage2/home/wangchenfei/annotations/mm9/mm9.len $1.Q10.nucleo.73bp.bw


## Fragment sampling
head -1000000 /mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/data.raw/nucleosome/TFKD/$1/$1.Q10.uniq.bed > $1.Q10.frag.bed

## Rotational positioning
python ../nucpos.py /mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/data.raw/nucleosome/TFKD/$1/$1.Q10.uniq.bed # Rotational positioning

## Statistical positionng
python ../nucpos.py /mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/data.raw/nucleosome/TFKD/$1/$1.Q10.uniq.bed # Statistical positioning

## Genome covergage
intersectBed -wa -a /mnt/Storage2/home/wangchenfei/annotations/mm9/mm9.200bp.bed -b /mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/data.processed/nucleosome/individual/$1/$1.Q10.uniq.nucleo.73bp.bed -u | wc -l
