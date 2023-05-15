#! /usr/bin/bash

# 1 构建STAR比对的索引文件，已经有的不需要再创建，需要注意文库长度是否相同
core=10
output_folder="/root/wangje/Reference/Star_index/Human/GeneCode/"
fa="/root/wangje/Reference/Homo_sapiens/GeneCode/hg38/Sequence/GRCh38.p13.genome.fa"
gtf="/root/wangje/Reference/Homo_sapiens/GeneCode/hg38/Annotation/Genes/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"


STAR \
    --runMode genomeGenerate \
    --runThreadN ${core} \
    --genomeDir ${output_folder} \
    --genomeFastaFiles  ${fa}\
    --sjdbGTFfile ${gtf} \
    --sjdbOverhang 149 #读段长度: 后续回帖读段的长度, 如果读长是PE 100， 则该值设为100-1=99

# 参数：
    # --runMode genomeGenerate \ #让STAR执行基因组索引的生成工作
    # --runThreadN 10 \ #构建运行使用的线程数
    # --genomeDir  #构建好的参考基因组存放的位置，最好是单独建立的一个文件夹
    # --genomeFastaFiles  #fasta文件（参考基因组序列文件）
    # --sjdbGTFfile  # gtf文件（基因注释文件）
    # --sjdbOverhang 149 #读段长度: 后续回帖读段的长度, 如果读长是PE 100， 则该值设为100-1=99

# 2 进行比对
STAR \
    --runThreadN $core \
    --outSAMstrandField intronMotif\
    --genomeDir $index --readFilesCommand zcat \
    --readFilesIn $input_folder/$sample*gz\
    --outFileNamePrefix $STAR_output_folder/$sample \
    --genomeLoad LoadAndKeep

# 参数
