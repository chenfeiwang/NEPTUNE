##################### Provide genome annotation files before you run this script #####################
# BWA reference
bwa_ref="/mnt/Storage/data/BWA/mm9.fa"
# chromosome length file
chr_bed="/mnt/Storage2/home/wangchenfei/annotations/mm9/chromInfo_mm9.txt"
# 200bp bins to calculate genome coverage
bin_bed="/mnt/Storage2/home/wangchenfei/annotations/mm9/mm9.200bp.bed"
# input PE fastqs
R1_fq=$1.1.fq
R2_fq=$1.2.fq

#####################  Step1: Quality control, for single sample running #####################
## mapping PE-seq MNase-seq data and filtering reads
bwa mem -t 10 $bwa_ref $R1_fq $R2_fq > $1.sam
samtools view -bht $chr_bed $1.sam > $1.bam
bamtools filter -mapQuality ">=10" -in $1.bam -out $1.Q10.bam

## uniq the mapped nucleosome reads
python NEPTUNE/NEPTUNE_quality_control.py uniq $1.Q10.bam mm9

## Tag extension and pileup, 37bp
python NEPTUNE/NEPTUNE_quality_control.py extend $1.Q10.uniq.bed 37
sort -k1,1 -k2,2n $1.Q10.uniq.nucleo.37bp.bed > $1.Q10.uniq.nucleo.37bp.tmp
mv $1.Q10.uniq.nucleo.37bp.tmp $1.Q10.uniq.nucleo.37bp.bed
genomeCoverageBed -i $1.Q10.uniq.nucleo.37bp.bed -g $chr_bed -bga -scale $2 > $1.Q10.nucleo.37bp.bdg
bedGraphToBigWig $1.Q10.nucleo.37bp.bdg $chr_bed $1.Q10.nucleo.37bp.bw

## Tag extension and pileup, 73bp
python NEPTUNE/NEPTUNE_quality_control.py extend $1.Q10.uniq.bed 73
sort -k1,1 -k2,2n $1.Q10.uniq.nucleo.73bp.bed > $1.Q10.uniq.nucleo.73bp.tmp
mv $1.Q10.uniq.nucleo.73bp.tmp $1.Q10.uniq.nucleo.73bp.bed
genomeCoverageBed -i $1.Q10.uniq.nucleo.73bp.bed -g $chr_bed -bga -scale $2 > $1.Q10.nucleo.73bp.bdg
bedGraphToBigWig $1.Q10.nucleo.73bp.bdg $chr_bed $1.Q10.nucleo.73bp.bw

## Fragment sampling
head -1000000 $1.Q10.uniq.bed > $1.Q10.frag.bed

## Rotational positioning
python NEPTUNE/NEPTUNE_quality_control.py rotational $1.Q10.uniq.bed $bwa_ref # Rotational positioning

## Statistical positionng
python NEPTUNE/NEPTUNE_quality_control.py statistical $1.Q10.uniq.bed $bwa_ref # Statistical positioning

## Genome covergage
intersectBed -wa -a $bin_bed -b $1.Q10.uniq.nucleo.73bp.bed -u | wc -l


#####################  Step2: NDR and POS score calculation, for multiple sample calculation #####################
heatmapr -b data/mm9.promoter.bed -w MalePSperm.Q10.uniq.nucleo.37bp.bw,MaleP0.5.Q10.uniq.nucleo.37bp.bw,MaleP1.Q10.uniq.nucleo.37bp.bw,MaleP1.5.Q10.uniq.nucleo.37bp.bw,MaleP2.Q10.uniq.nucleo.37bp.bw,MaleP3.Q10.uniq.nucleo.37bp.bw,MaleP4.Q10.uniq.nucleo.37bp.bw,MaleP6.Q10.uniq.nucleo.37bp.bw,MaleP8.Q10.uniq.nucleo.37bp.bw,MaleP12.Q10.uniq.nucleo.37bp.bw, --name=MalePN_promoter --method=mean --s_wigindex=10 --dir --upstream=1000 --downstream=1000
heatmapr -b data/mm9.promoter.bed -w FemalePMII.Q10.uniq.nucleo.37bp.bw,FemaleP0.5.Q10.uniq.nucleo.37bp.bw,FemaleP1.Q10.uniq.nucleo.37bp.bw,FemaleP1.5.Q10.uniq.nucleo.37bp.bw,FemaleP2.Q10.uniq.nucleo.37bp.bw,FemaleP3.Q10.uniq.nucleo.37bp.bw,FemaleP4.Q10.uniq.nucleo.37bp.bw,FemaleP6.Q10.uniq.nucleo.37bp.bw,FemaleP8.Q10.uniq.nucleo.37bp.bw,FemaleP12.Q10.uniq.nucleo.37bp.bw --name=FemalePN_promoter --method=mean --s_wigindex=10 --dir --upstream=1000 --downstream=1000
python2.7 NEPTUNE/NEPTUNE_profile_score.py MalePN 10
python2.7 NEPTUNE/NEPTUNE_profile_score.py FemalePN 10

python2.7 NEPTUNE/NEPTUNE_profile_ZGA.py MalePN 10
python2.7 NEPTUNE/NEPTUNE_profile_ZGA.py FemalePN 10
















