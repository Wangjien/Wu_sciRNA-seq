#! /usr/bin/bash

# 循环调用文件中的05_removeDuplicateUmi.py
input_folder="root/wangje/Project/Yin/LB230410/new_output/o4_splitsam"
sample_prefix="LB230410-"
output_folder="root/wangje/Project/Yin/LB230410/new_output/05"

for i in {1..96}
do 
    python /root/wangje/Project/Yin/LB230410/new_output/05.py ${input_folder} ${sample_prefix}_${i}_filterAndSort.sam ${output_folder}
done

