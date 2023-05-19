import pandas as pd
import numpy as np
import os
import gzip
import sys
import time
import pyfastx

# 读入文件
def Extract_01(input_folder, output_folder, sample_id, ligation_file, RT_file):
    Read1 = os.path.join(input_folder,f'{sample_id}_1.fastq.gz')
    Read2 = os.path.join(input_folder,f'{sample_id}_2.fastq.gz')
    output_file = os.path.join(output_folder,f'{sample_id}_barcode.fastq.gz')
    f3 = gzip.open(output_file, 'wb')
    
    ligation_list =  pd.read_csv(ligation_file, sep = '\t', header=None, names=['ligation'])
    ligation_list = ligation_list["ligation"].tolist()
    
    RT_list = pd.read_csv(RT_file, header=None, names=['RT'])
    RT_list = RT_list["RT"].tolist()
     
    names = []
    first_line = []
    for name,seq,qual in pyfastx.Fastq(Read1, build_index=False):
        tmp_lig = seq[0:10]
        if tmp_lig in ligation_list:
            tmp_RT = seq[len(tmp_lig) + 14 : len(tmp_lig) + 24]
            if tmp_RT in RT_list:
                umi = seq[len(tmp_lig)+6: len(tmp_lig) + 14]
                names.append(name) 
                first_line.append("@"+ name + ',' + tmp_lig + tmp_RT + "," + umi)
    names_dict = dict(zip(names, first_line))
    
    count = 0
    for name,seq,qual in pyfastx.Fastq(Read2, build_index=False):
        if name in names:
            f3.write((names_dict[name] +"\n").encode())
            f3.write((seq + "\n").encode())
            f3.write(("+" + "\n").encode())
            f3.write((qual + "\n").encode())
            count +=1
            percent = (count/len(names_dict))*100
            print('\r' + "#" * round(percent) + ">>> %d%% (%d/%d)" %(percent,count,len(names_dict)),end="")
    f3.close()    
        
if __name__ == '__main__':
    
    input_folder="/root/wangje/Project/吴霞/SRR7827964"
    output_folder = "/root/wangje/Project/吴霞/SRR7827964"
    sample_id = "SRR7827964"
    ligation_file = "/root/wangje/Project/吴霞/SRR7827964/ligation_file.txt"
    RT_file = "/root/wangje/Project/吴霞/SRR7827964/RT_file.txt"
    Extract_01(
        input_folder=input_folder,
        output_folder= output_folder,
        sample_id=sample_id,
        ligation_file=ligation_file,
        RT_file=RT_file
    )
