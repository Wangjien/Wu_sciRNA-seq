#!/usr/bin/bash

input_folder="/root/wangje/Project/Yin/LB230410/new_output/03_STAR_output"
output_folder="/root/wangje/Project/Yin/LB230410/new_output/04_splitsam"
for i in {1..96}
do
    echo ">>>>>>>>>> sample: BC230502-${i}"
    samtools view -bh -q 30 -F 4 $input_folder/LB230410-${i}-Aligned.out.sam|samtools sort -@ 10 -|samtools view -h ->$output_folder/LB230410-${i}_filterAndSort.sam
done

# 参数：
# -b /-h :-b indicates that the output should be in the BAM format (binary SAM), while -h indicates the output should be in the SAM format.
# -F 4：统计map 上的 reads总数
# -q This parameter is used to filter alignments based on their mapping quality. Alignments with mapping quality below the specified threshold will be excluded. 