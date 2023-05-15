import os
import sys
import time
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from  Levenshtein import distance

"""

"""
def rm_dup_samfile(input_folder, sample_prefix,output_folder):
    input_file = os.path.join(input_folder, f"{sample_prefix}.sam")
    output_file = os.path.join(output_folder, f"{sample_prefix}_rmDup.sam")
    f2 = open(output_file, 'w')
    all_counts = []
    with open(input_file) as f1:
        for l1 in f1.readlines():
            all_counts.append(l1)
    all_num = len(all_counts)
    all_counts = None
        
    pre_barcode = []
    pre_line = []
    unique_id = []
    pre_chrom = 0
    pre_site = 0
    dup = False
    
    pre_dup_num = 0
    cur_dup_num = 0
    count = 0
    
    
    for line in open(input_file):
        count +=1
        if line.startswith("@"):
            # print(line.strip())
            f2.write(line)
            continue
        else:
            name = line.strip().split("\t")[0].split(",")
            barcode_UMI = name[1]+name[2]
            chrom_num = line.strip().split("\t")[2]
            start_site = line.strip().split("\t")[3]
            
            if((start_site == pre_site) and (chrom_num == pre_chrom)):
                dup = False
                for each_barcode in pre_barcode:
                    if each_barcode == barcode_UMI:
                        dup = True
                        break
                    if dup == False:
                        pre_dup_num = cur_dup_num
                        cur_dup_num = 1
                        f2.write(line)
                        # print(line)
                    else:
                        cur_dup_num += 1
            else:
                pre_dup_num = cur_dup_num
                cur_dup_num = 1
                f2.write(line)
                # print(line)
                pre_chrom = chrom_num
                pre_site = start_site
                pre_barcode = set()
                pre_barcode.add(barcode_UMI)
        percent = (int(count)/int(all_num)) * 100
        # time.sleep(0.0001)
        print("\r" + "#"* int(percent),"->[%s] %d%% (%d/%d)" %(sample_prefix,percent,count,all_num), end = "")
    f2.close()
    
    
if __name__ == '__main__':
    input_folder = sys.argv[1]
    sample_prefix = sys.argv[2]
    output_folder = sys.argv[3]
    
    # input_folder = "/root/wangje/Project/吴霞/Data/04_filterAndSort/"
    # sample_prefix = "BC230502-10_filterAndSort"
    # output_folder = "/root/wangje/Project/吴霞/Data/05_removeDDuplicateUMI/"
    rm_dup_samfile(
        input_folder=input_folder,
        sample_prefix=sample_prefix,
        output_folder=output_folder
    )
    
      
    

                
                    
            
            
            
        
    
    


