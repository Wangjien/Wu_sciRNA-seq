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
## 的那个文件比对
# STAR \
#     --runThreadN $core \
#     --outSAMstrandField intronMotif\
#     --genomeDir $index --readFilesCommand zcat \
#     --readFilesIn $input_folder/$sample*gz\
#     --outFileNamePrefix $STAR_output_folder/$sample \
#     --genomeLoad LoadAndKeep

# 参数
# --outSAMstrandField  非链特异性RNA-seq，同时为了保证能和Cufflinks兼容
# --genomeDir 索引文件存放路径 : 必须先创建文件夹
# --readFilesCommand 如果输入格式是gz结尾，那么需要加上zcat, 否则会报错
# --readFilesIn 输入文件
# --outFileNamePrefix 输出文件qianz
# --genomeLoad 多个比对共享内存中的基因组索引 减小内存的使用量

## 并发，加快运行速度，注意内存占用
core=4
index="/root/wangje/Reference/Star_index/Human/GeneCode/"
STAR_output_folder="/root/wangje/Project/吴霞/Data/03_STAR_output"
input_folder="/root/wangje/Project/吴霞/Data/02_trim"

tempfifo="my_temp_fifo"
mkfifo ${tempfifo}
rm -rf ${tempfifo}
exec ${core}<>${tempfifo}

for i in {1..96}
do
{
    STAR \
    --runThreadN $core \
    --outSAMstrandField intronMotif\
    --genomeDir $index --readFilesCommand zcat \
    --readFilesIn $input_folder/BC230502-${i}_R2_barcode_trimmed.fq.gz\
    --outFileNamePrefix $STAR_output_folder/BC230502-${i} \
    --genomeLoad LoadAndKeep
}& 
done  >&4 
wait
exec 4>&-
