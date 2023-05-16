#! /usr/bin/bash
echo "提取Read1中的序列,添加到满足条件的Read2中并输出"
# input_folder="/root/wangje/Project/吴霞" 
# output_folder="/root/wangje/Project/吴霞/Data" 
# sample_id="BC230502-1"
# ligation_barcode_file="/root/wangje/Project/吴霞/script/ligation_barcode.txt"
# RT_Barcode_file="/root/wangje/Project/吴霞/script/RT_barcode.txt"
# p7_file="/root/wangje/Project/吴霞/script/sample_p7.txt"

# tmpfifo='my_tmp_fifo'
# mkfifo ${tmpfifo}
# exec 8 <> ${tmpfifo}
# rm -f ${tmpfifo}

# for i in {1..96}
#     do {
#         python /root/wangje/Project/吴霞/script/01_根据umi提取符合条件的barcode并输出文件.py ${input_folder}/$sample_id $output_folder $sample_id $ligation_barcode_file $RT_Barcode_file $p7_file

#     }&
# done >& 8
# wait
# exec 8>&-     

#################################################
# 从Read1中提取符合条件的序列并且输出文件
#################################################

input_folder="/root/wangje/Project/吴霞" 
output_folder="/root/wangje/Project/吴霞/Data" 
ligation_barcode_file="/root/wangje/Project/吴霞/script/ligation_barcode.txt"
RT_Barcode_file="/root/wangje/Project/吴霞/script/RT_barcode.txt"
p7_file="/root/wangje/Project/吴霞/script/sample_p7.txt"
for i in {1..96}
do 
    sample_id="BC230502-${i}"
    echo "$sample_id"
    python /root/wangje/Project/吴霞/script/01_根据umi提取符合条件的barcode并输出文件.py ${input_folder}/$sample_id $output_folder $sample_id $ligation_barcode_file $RT_Barcode_file $p7_file
done

####################################################
# 将sam文件按照barcode类型进行分割
####################################################
echo ">>>>>>>>>> split samfile by barcode <<<<<<<<<<"
sam_file="/home/data/vip9t22/Project/BC2305/BC230502-" 
barcode_file="/home/data/vip9t22/Project/BC2305/combined_barcode.txt" 
output_folder="/home/data/vip9t22/Project/BC2305/BC230502-"
cutoff=20
for i in {1..96}
do 
    python 06.py ${sam_file}${i}_filterAndSort_rmDup.sam ${barcode_file} ${output_folder}${i} $cutoff
done

###################################################
# gene 计数
###################################################
gtf_file="/root/wangje/Reference/Homo_sapiens/GeneCode/hg38/Annotation/Genes/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"
input_folder="/root/wangje/Project/吴霞/Data/06_splitSAM"
sample_ID="/root/wangje/Project/吴霞/Data/06_splitSAM"
core=10

for i in {1..96}
do 
    python  /root/wangje/Project/吴霞/script/07.py ${gtf_file} ${input_folder}/BC230502-${i} ${sample_ID}/BC230502-${i}/BC230502-${i}_filterAndSort_rmDup.sample_list.txt ${core}
done



