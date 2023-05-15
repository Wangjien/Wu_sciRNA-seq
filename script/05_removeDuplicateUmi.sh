#! /usr/bin/bash

# 循环调用文件中的05_removeDuplicateUmi.py
input_folder="/root/wangje/Project/吴霞/Data/04_filterAndSort",
sample_prefix="BC230502-",
output_folder="/root/wangje/Project/吴霞/Data/05_removeDDuplicateUMI"

for i in {1..96}
do 
    python /root/wangje/Project/吴霞/script/05_removeDuplicateUmi.py ${input_folder} ${sample_prefix}${i}_filterAndSort.sam ${output_folder}
done

