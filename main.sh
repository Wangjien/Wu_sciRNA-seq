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

input_folder="/root/wangje/Project/Yin/LB230410/Data" 
output_folder="/root/wangje/Project/Yin/LB230410/new_output" 
ligation_barcode_file="/root/wangje/Project/吴霞/script/ligation_barcode.txt"
RT_Barcode_file="/root/wangje/Project/吴霞/script/RT_barcode.txt"
p7_file="/root/wangje/Project/吴霞/script/sample_p7.txt"
for i in {1..96}
do 
    sample_id="LB230410-${i}_LB230410-${i}"
    echo "$sample_id"
    python /root/wangje/Project/Yin/LB230410/01.py ${input_folder} $output_folder $sample_id $ligation_barcode_file $RT_Barcode_file $p7_file
done

####################################################
# 将sam文件按照barcode类型进行分割
####################################################
echo ">>>>>>>>>> split samfile by barcode <<<<<<<<<<"
# 设置最大文件输出数目，linux命令
ulimit -n 100000
# sam_file="/home/data/vip9t22/Project/BC2305/BC230502-" 
# barcode_file="/home/data/vip9t22/Project/BC2305/combined_barcode.txt" 
# output_folder="/home/data/vip9t22/Project/BC2305/BC230502-"

sam_file="/root/wangje/Project/吴霞/Data/05_removeDDuplicateUMI/BC230502-"
barcode_file="/root/wangje/Project/吴霞/Data/combined_barcode.txt"
output_folder="/root/wangje/Project/吴霞/Data/06_splitSam_new/BC230502-"
cutoff=3
for i in {1..96}
do 
    python /root/wangje/Project/吴霞/script/06_splitSAM.py ${sam_file}${i}_filterAndSort_rmDup.sam ${barcode_file} ${output_folder}${i} $cutoff
done

###################################################
# gene 计数
###################################################
gtf_file="/root/wangje/Reference/Homo_sapiens/GeneCode/hg38/Annotation/Genes/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"
# input_folder="/root/wangje/Project/吴霞/Data/06_splitSAM"
input_folder="/root/wangje/Project/吴霞/Data/06_splitSam_new"
# sample_ID="/root/wangje/Project/吴霞/Data/06_splitSAM"
sample_ID="/root/wangje/Project/吴霞/Data/06_splitSam_new"
core=40

for i in {1..96}
do 
    python  /root/wangje/Project/吴霞/script/07.py ${gtf_file} ${input_folder}/BC230502-${i} ${sample_ID}/BC230502-${i}/BC230502-${i}_filterAndSort_rmDup.sample_list.txt ${core}
done

## 并发，提高运行速度（有问题需改进）

# gtf_file="/root/wangje/Reference/Homo_sapiens/GeneCode/hg38/Annotation/Genes/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz"
# input_folder="/root/wangje/Project/吴霞/Data/06_splitSAM"
# sample_ID="/root/wangje/Project/吴霞/Data/06_splitSAM"
# core=10

# tempfifo="my_temp_fifo"
# mkfifo ${tempfifo}
# rm -rf ${tempfifo}
# exec 4<>${tempfifo}

# for i in {1..96}
# do 
# {
#     python  /root/wangje/Project/吴霞/script/07.py ${gtf_file} ${input_folder}/BC230502-${i} ${sample_ID}/BC230502-${i}/BC230502-${i}_filterAndSort_rmDup.sample_list.txt ${core}
# }&
# done >&4 
# wait
# exec 4>&-

#######################################################
# 转化第7步的结果，转换为Seurat Read10x函数的输入文件
#######################################################
pwd="/root/wangje/Project/吴霞/Data/06_splitSAM"
script="/root/wangje/Project/吴霞/script/08_转换第7步的结果.py"

for i in {1..96}
do  
    echo ${pwd}/BC230502-${i}
    python $script $pwd/BC230502-${i}
done

#########################################################
# 合并文件
#########################################################
for i in {2..96}
do 
    # cd /root/wangje/Project/吴霞/Data/06_splitSAM/BC230502-${i}
    cat /root/wangje/Project/吴霞/Data/06_splitSAM/BC230502-${i}/*count > /root/wangje/Project/吴霞/Data/06_splitSAM/BC230502-${i}/outfile/count.MM
    cat /root/wangje/Project/吴霞/Data/06_splitSAM/BC230502-${i}/*report > /root/wangje/Project/吴霞/Data/06_splitSAM/BC230502-${i}/outfile/report.MM
    cp /root/wangje/Project/吴霞/Data/06_splitSAM/BC230502-${i}/*txt  /root/wangje/Project/吴霞/Data/06_splitSAM/BC230502-${i}/outfile/
done 



for i in {2..96}
do  
    cat /root/wangje/Project/Yin/LB230410/new_output/06/test${i}/*.count > /root/wangje/Project/Yin/LB230410/new_output/06/test${i}/outfile/count.MM
    cat /root/wangje/Project/Yin/LB230410/new_output/06/test${i}/*.report > /root/wangje/Project/Yin/LB230410/new_output/06/test${i}/outfile/report.MM
    cp /root/wangje/Project/Yin/LB230410/new_output/06/test${i}/*.txt  /root/wangje/Project/Yin/LB230410/new_output/06/test${i}/outfile/
done













