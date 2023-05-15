#!/usr/bin/bash

input_folder="/root/wangje/Project/吴霞/Data/03_STAR_output"
output_folder="/root/wangje/Project/吴霞/Data/04_filterAndSort"
for i in {1..96}
do
    echo ">>>>>>>>>> sample: BC230502-${i}"
    samtools view -bh -q 30 -F 4 $input_folder/BC230502-${i}Aligned.out.sam|samtools sort -@ 10 -|samtools view -h ->$output_folder/BC230502-${i}_filterAndSort.sam
done

# 参数：
# -b /-h :-b indicates that the output should be in the BAM format (binary SAM), while -h indicates the output should be in the SAM format.
# -F 4：统计map 上的 reads总数
# -q This parameter is used to filter alignments based on their mapping quality. Alignments with mapping quality below the specified threshold will be excluded.

