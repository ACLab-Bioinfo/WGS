# ---------------------------->>>>>>>>>>
# GWS
# ---------------------------->>>>>>>>>>

# pipeline
# 1. QC
# 原始测序数据 ——> QC/过滤低质量read数据
# 2. Data Preprocess
# reads Alignment ——> sort ——> redup ——> 局部重比对 ——> 碱基质量重校正(BQSR) 
# 3. 变异检测
# 单样本：变异质控和过滤(VQSR)(VCF)
# 多样本：Merge(gVCF)--Joint Genotype(gVCF)--变异质控和过滤(VQSR)(VCF)


#!/bin/bash
## 血的教训，切勿直接套用别人写好的代码，不同版本用法可能区别很大
## 本流程适用于 The Genome Analysis Toolkit (GATK) v4.1.2.0
## 这个流程假设你只有一个样本，这个样本只有一对用Illumina测序仪测序的PE fastq数据文件。
# qsub run.lab5.wgs_single.sh
# bash lab5.wgs_single.sh lab4_1.fastq.gz lab4_2.fastq.gz SRR1264357 WXS.chr22 HCT-116 IDX1 HCT116.WXS.chr22
# bash lab5.wgs_single.sh WXS-Case_R1.fq WXS-Case_R2.fq Case1 WXS Case1 IDX1 Case1.WXS
hostname
date
echo "** Starting to call SNP/indel from fastq for a single sample **"

# 由conda所装
# which trimmomatic 
# /public/workspace/Course/NGS/anaconda/bin/trimmomatic
# which gatk
# /public/workspace/Course/NGS/anaconda/bin/gatk

# 通过module模块调入bwa和samtools，随后的使用无需制定路径
# module load bwa.kit_0.7.15 会自动加载bwa 0.7.15和samtools 1.3
# bytlib load bwa.kit_0.7.15
# bytlib load samtools-1.9
# 一些软件和工具的路径, 根据实际
# bwa=/bio-apps/rhel7/bwa.kit_0.7.15/bwa
# samtools=/bio-apps/rhel7/bwa.kit_0.7.15/samtools
# samtools=/bio-apps/rhel7/samtools-1.9/bin/samtools
# trimmomatic=/your_path_to/Trimmomatic/0.36/trimmomatic-0.36.jar
# gatk=/your_path_to/gatk/4.0.3.0/gatk

# reference
# reference=/public/workspace/Course/NGS/Reference/BWA_Index/Homo_sapiens.GRCh38.dna.chromosome.22.fa
reference=/public/workspace/shaojf/Course/NGS/Reference/BWA_Index/Homo_sapiens.GRCh38.dna.primary_assembly.fa
GATK_bundle=/public/workspace/shaojf/Course/NGS/Reference/ftp.broadinstitute.org/bundle/hg38

## 这一步不是必须的，取决于GATK_bundle中的这4份文件是否已经有建索引，如没有再执行
#$gatk IndexFeatureFile --feature-file $GATK_bundle/hapmap_3.3.hg38.vcf  
#$gatk IndexFeatureFile --feature-file $GATK_bundle/1000G_omni2.5.hg38.vcf  
#$gatk IndexFeatureFile --feature-file $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf  
#$gatk IndexFeatureFile --feature-file $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/dbsnp_146.hg38.vcf  

## shell执行参数
fq1=$1
fq2=$2
RGID=$3  ## Read Group，一般用Lane ID代替
library=$4  ## 测序文库编号
sample=$5  ## 样本ID
pu=$6	## Read-Group platform unit
outdir=$7  ## 输出目录的路径

# fq1=WXS-Case_R1.fq
# fq2=WXS-Case_R2.fq
# RGID=Case1
# library=WXS
# sample=Case1
# pu=IDX1
# outdir=Case1.WXS

## 按样本设置目录
outdir=${outdir}/${sample}

## 通过fastq1获得fastq的前缀名字，这里假设了原始的fastq1和fastq2有相同的前缀名字
fq_file_name=`basename $fq1`
fq_file_name=${fq_file_name%%.*}

# output diretory
if [ ! -d $outdir/cleanfq ]
then mkdir -p $outdir/cleanfq
fi

if [ ! -d $outdir/bwa ]
then mkdir -p $outdir/bwa
fi

if [ ! -d $outdir/gatk ]
then mkdir -p $outdir/gatk
fi

## 使用Trimmomatic对原始数据进行质控，ILLUMINACLIP中的一个关键参数 keepBothReads设为True。
# time java -jar ${trimmomatic} PE \
time trimmomatic PE -threads 20 \
	$fq1 $fq2 \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz ${fq_file_name}.unpaired.1.fq.gz \
	$outdir/cleanfq/${fq_file_name}.paired.2.fq.gz ${fq_file_name}.unpaired.2.fq.gz \
	ILLUMINACLIP:/public/workspace/Course/NGS/Tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:True \
	SLIDINGWINDOW:5:15 LEADING:5 TRAILING:5 MINLEN:50 && echo "** fq QC done **"

## 使用bwa mem完成数据比对，bwa mem对任何长度大于40bp小于2000bp的read都是非常有效的; PL:ILLUMINA是我默认的
time bwa mem -t 24 -M -Y -R "@RG\tID:$RGID\tPL:ILLUMINA\tPU:$pu\tLB:$library\tSM:$sample" $reference \
	$outdir/cleanfq/${fq_file_name}.paired.1.fq.gz $outdir/cleanfq/${fq_file_name}.paired.2.fq.gz | \
	samtools view -Sb - > $outdir/bwa/${sample}.bam
echo "** BWA MEM done **"
# BWA比对后输出的BAM文件没有顺序
# -@，用于设定排序时的线程数
# -m，限制排序时最大的内存消耗
# -O 指定输出为bam格式
# -o 是输出文件的名字
time samtools sort -@ 8 -m 4G -O bam -o $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.bam
echo "** sorted raw bamfile done **"

## 这一步不是必须的 
# time samtools index $outdir/bwa/${sample}.sorted.bam && echo "** ${sample}.sorted.bam index done **"
# picard AddOrReplaceReadGroups \
# 	I=$outdir/bwa/${sample}.sorted.bam \
# 	O=$outdir/bwa/${sample}.sorted.replaceRG.bam \
# 	RGID=$RGID RGLB=$library RGPL=Illumina RGPU=$pu RGSM=$sample

## 标记重复序列 
gatk MarkDuplicates \
	-I $outdir/bwa/${sample}.sorted.bam \
	-M $outdir/bwa/${sample}.markdup_metrics.txt \
	-O $outdir/bwa/${sample}.sorted.markdup.bam
echo "** ${sample}.sorted.bam MarkDuplicates done **"
  
## 为${sample}.sorted.markdup.bam构建Index，这是继续后续步骤所必须的
time samtools index $outdir/bwa/${sample}.sorted.markdup.bam
echo "** ${sample}.sorted.markdup.bam index done **"

## 多样本需要进行局部重比对（GATK）
# 第一步，RealignerTargetCreator ，目的是定位出所有需要进行序列重比对的目标区域
# 第二步，IndelRealigner，对所有在第一步中找到的目标区域运用算法进行序列重比对，最后得到捋顺了的新结果
##

## 执行BQSR
# BQSR（Base Quality Score Recalibration）主要是通过机器学习的方法构建测序碱基的错误率模型，
# 然后对这些碱基的质量值进行相应的调整。
# 第一步，BaseRecalibrator，这里计算出了所有需要进行重校正的read和特征值，然后把这些信息输出为一份校准表文件（sample_name.recal_data.table）
# 第二步，PrintReads，这一步利用第一步得到的校准表文件（sample_name.recal_data.table）重新调整原来BAM文件中的碱基质量值，并使用这个新的质量值重新输出一份新的BAM文件。

## [注]Does your vcf file have an index? GATK4 does not support on the fly indexing of VCFs anymore.
time gatk BaseRecalibrator \
	-R $reference \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	--known-sites $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
	--known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
	--known-sites $GATK_bundle/dbsnp_146.hg38.vcf \
	-O $outdir/bwa/${sample}.sorted.markdup.recal_data.table
echo "** ${sample}.sorted.markdup.recal_data.table done **" 

time gatk ApplyBQSR \
	--bqsr-recal-file $outdir/bwa/${sample}.sorted.markdup.recal_data.table \
	-R $reference \
	-I $outdir/bwa/${sample}.sorted.markdup.bam \
	-O $outdir/bwa/${sample}.sorted.markdup.BQSR.bam
echo "** ApplyBQSR done **"

## 为${sample}.sorted.markdup.BQSR.bam构建Index，这是继续后续步骤所必须的
time samtools index $outdir/bwa/${sample}.sorted.markdup.BQSR.bam
echo "** ${sample}.sorted.markdup.BQSR.bam index done **"


## 对于单个样本来说，有四个完成变异检测的方式，结果是一样的，可以按照需要挑选一种(以下默认第一种)。
## 第一，直接调用HaplotypeCaller输出样本VCF，面对较大的输入文件时，速度较慢
# HaplotypeCaller和那些直接应用贝叶斯推断的算法有所不同，
# 它会先推断群体的单倍体组合情况，计算各个组合的几率，
# 然后根据这些信息再反推每个样本的基因型组合。
# 因此它不但特别适合应用到群体的变异检测中，
# 而且还能够依据群体的信息更好地计算每个个体的变异数据和它们的基因型组合。

time gatk HaplotypeCaller \
	-R $reference \
	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
	-O $outdir/gatk/${sample}.HC.vcf.gz
echo echo "** ${sample}.HC.vcf.gz done ** "

# ## 第二，输出这个样本每个染色体的vcf，然后在合并所有的染色体结果，目的是提高速度，这不是必须的，仅是通过分染色体获得速度的提升
# chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM )
# for i in ${chrom[@]}; do
# 	time $gatk HaplotypeCaller \
# 	  -R $reference/Homo_sapiens_assembly38.fasta \
#       -I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#       -L $i \
#       -O $outdir/gatk/${sample}.HC.${i}.vcf.gz && echo "** ${sample}.HC.${i}.vcf.gz done ** " &
# done && wait
# merge_vcfs=""
# for i in ${chrom[@]}; do
#     merge_vcfs=${merge_vcfs}" -I $outdir/gatk/${sample}.HC.${i}.vcf.gz \\"\n
# done && time $gatk MergeVcfs ${merge_vcfs} -O $outdir/gatk/${sample}.HC.vcf.gz && echo "** MergeVcfs done **"

# ## 第三，先输出样本的全gVCF，再进行GenotypeGVCFs，这个方式在单样本情况下不是必须的，但是多样本的标配，面对较大的输入文件时，速度较慢
# time $gatk HaplotypeCaller \
#     --emit-ref-confidence GVCF \
# 	-R $reference/Homo_sapiens_assembly38.fasta \
# 	-I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
# 	-O $outdir/gatk/${sample}.HC.g.vcf.gz && echo "** GVCF ${sample}.HC.g.vcf.gz done **"
# time $gatk GenotypeGVCFs \
# 	-R $reference/Homo_sapiens_assembly38.fasta \
# 	-V $outdir/gatk/${sample}.HC.g.vcf.gz \
# 	-O $outdir/gatk/${sample}.HC.vcf.gz && echo "** ${sample}.HC.vcf.gz done ** "

# ## 第四，输出每个染色体的gvcf，然后对每个染色体单独进行GenotypeGVCFs，目的是提高速度，这不是必须的，仅是通过分染色体获得速度的提升
# chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM )
# for i in ${chrom[@]}; do
# 	time $gatk HaplotypeCaller \
#       --emit-ref-confidence GVCF \
# 	    -R $reference/Homo_sapiens_assembly38.fasta \
#       -I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam \
#       -L $i \
# 		-O $outdir/gatk/${sample}.HC.${i}.g.vcf.gz && \
# 	time $gatk GenotypeGVCFs \
# 		-R $reference/Homo_sapiens_assembly38.fasta \
# 		-V $outdir/gatk/${sample}.HC.${i}.g.vcf.gz \
# 		-O $outdir/gatk/${sample}.HC.${i}.vcf.gz && echo "** ${sample}.HC.${i}.vcf.gz done ** " & 	
# done && wait
# merge_vcfs=""
# for i in ${chrom[@]}; do
#     merge_vcfs=${merge_vcfs}" -I $outdir/gatk/${sample}.HC.${i}.vcf.gz \\"\n
# done && time $gatk MergeVcfs ${merge_vcfs} -O $outdir/gatk/${sample}.HC.vcf.gz && echo "** MergeVcfs done **"

## VQSR, 由于评价SNP和Indel质量高低的标准是不同的，因此，需要分SNP和Indel这两种不同的模式，分别进行质控
## 首先是SNP mode
# --max-gaussians 4 默认是8，往下调有利于小数据
# --rscript-file $outdir/gatk/${sample}.HC.snps.plots.R \
time gatk VariantRecalibrator \
	-R $reference \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf \
	-resource:omini,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATK_bundle/dbsnp_146.hg38.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode SNP \
	--max-gaussians 4 \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
	-O $outdir/gatk/${sample}.HC.snps.recal
time gatk ApplyVQSR \
	-R $reference \
	-V $outdir/gatk/${sample}.HC.vcf.gz \
	--tranches-file $outdir/gatk/${sample}.HC.snps.tranches \
	--recal-file $outdir/gatk/${sample}.HC.snps.recal \
	--truth-sensitivity-filter-level 99.0 \
	-mode SNP \
	-O $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz
echo "** SNPs VQSR done **"
# [June 4, 2019 6:28:34 PM CST] org.broadinstitute.hellbender.tools.walkers.vqsr.VariantRecalibrator done. Elapsed time: 15.45 minutes.
# Runtime.totalMemory()=2023227392
# Tool returned:
# true
# Exception in thread "Thread-1" htsjdk.samtools.util.RuntimeIOException: java.nio.file.NoSuchFileException: /tmp/Rlib.6794369076579920831
#         at htsjdk.samtools.util.IOUtil.recursiveDelete(IOUtil.java:1346)
#         at org.broadinstitute.hellbender.utils.io.IOUtils.deleteRecursively(IOUtils.java:1061)
#         at org.broadinstitute.hellbender.utils.io.DeleteRecursivelyOnExitPathHook.runHooks(DeleteRecursivelyOnExitPathHook.java:56)
#         at java.lang.Thread.run(Thread.java:745)
# Caused by: java.nio.file.NoSuchFileException: /tmp/Rlib.6794369076579920831
#         at sun.nio.fs.UnixException.translateToIOException(UnixException.java:86)
#         at sun.nio.fs.UnixException.rethrowAsIOException(UnixException.java:102)
#         at sun.nio.fs.UnixException.rethrowAsIOException(UnixException.java:107)
#         at sun.nio.fs.UnixFileAttributeViews$Basic.readAttributes(UnixFileAttributeViews.java:55)
#         at sun.nio.fs.UnixFileSystemProvider.readAttributes(UnixFileSystemProvider.java:144)
#         at sun.nio.fs.LinuxFileSystemProvider.readAttributes(LinuxFileSystemProvider.java:99)
#         at java.nio.file.Files.readAttributes(Files.java:1737)
#         at java.nio.file.FileTreeWalker.getAttributes(FileTreeWalker.java:219)
#         at java.nio.file.FileTreeWalker.visit(FileTreeWalker.java:276)
#         at java.nio.file.FileTreeWalker.walk(FileTreeWalker.java:322)
#         at java.nio.file.Files.walkFileTree(Files.java:2662)
#         at java.nio.file.Files.walkFileTree(Files.java:2742)
#         at htsjdk.samtools.util.IOUtil.recursiveDelete(IOUtil.java:1344)
#         ... 3 more


## 然后是Indel mode
# --rscript-file $outdir/gatk/${sample}.HC.snps.indels.plots.R \
time gatk VariantRecalibrator \
	-R $reference \
	-V $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
	-resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode INDEL \
	--max-gaussians 2 \
	--tranches-file $outdir/gatk/${sample}.HC.snps.indels.tranches \
	-O $outdir/gatk/${sample}.HC.snps.indels.recal
time gatk ApplyVQSR \
	-R $reference \
	-V $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz \
	--truth-sensitivity-filter-level 99.0 \
	--tranches-file $outdir/gatk/${sample}.HC.snps.indels.tranches \
	--recal-file $outdir/gatk/${sample}.HC.snps.indels.recal \
	-mode INDEL \
	-O $outdir/gatk/${sample}.HC.VQSR.vcf.gz
echo "** SNPs and Indels VQSR (${sample}.HC.VQSR.vcf.gz finish) done **"

# 可以被删除清理的文件，这不是必须执行的
# 1) 对于比对文件只有最终的${sample}.sorted.markdup.BQSR.bam值得留下来
# 2) 对于VCF，可以保留.g.vcf.gz，原始HC.vcf.gz和完成质控的HC.VQSR.vcf.gz
# rm -f $outdir/bwa/${sample}.bam $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.sorted.markdup.bam*
# rm -f $outdir/gatk/${sample}.HC.snps.VQSR.vcf.gz rm -f $outdir/gatk/${sample}.HC.*.recal

echo "** finished calling SNP/indel from fastq for a single sample **"
date

#!/bin/bash
# run bash
bash lab5.wgs_single.sh lab4_1.fastq.gz lab4_2.fastq.gz SRR1264357 WXS.chr22 HCT-116 IDX1 HCT116.WXS.chr22
bash lab5.wgs_single.sh WXS-Case_R1.fq WXS-Case_R2.fq Case1 WXS Case1 IDX1 Case1.WXS
