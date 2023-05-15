import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

def samfile_barcode_count(samfile, barcode_file):
    """
    统计每种barcode测到的数目,单个文件和所有文件
    """
    # 生成最初的barcode计算问文件
    barcode_file = pd.read_csv(barcode_file, sep = '\t', header=None, names = ["barcode"])
    barcode_list = barcode_file["barcode"].tolist()
    barcode_conut = {k : 0 for k in barcode_list}
    
    # 统计文件中的barcode个数
    with open(samfile) as fw:
        for line in fw.readlines():
            if line.startswith("@"):
                continue
            else:
                line = line.strip().split("\t")
                barcode = line[0].split(",")[1]
                if barcode in barcode_conut:
                    barcode_conut[barcode] += 1
        return barcode_conut
def split_samfile(samfile, barcode_file, output_folder, cutoff):
    barcode_count = samfile_barcode_count(samfile, barcode_file)
    plot_name = samfile.split("_")[0]
    fig = plt.figure()
    plt.hist(barcode_count.values(), bins=100)
    plt.ylabel('frequency')
    plt.xlabel("Number of unique read")
    fig_output = os.path.join(output_folder,f"{plot_name}.png")
    
    fig.savefig(fig_output)
    
    # 写出每个文件的barcode分布
    read_dist = open(os.path.join(output_folder,f"{plot_name}_Barcode_distribution.txt"), 'w')
    for barcode, value in barcode_count.items():
        read_dist.write("%s \t %s" %(barcode, value))
    read_dist.close()
    
    # 根据cutoff对barcode进行过滤
    barcode_filtered = []
    for barcode in barcode_count:
        if barcode_count[barcode] >= cutoff:
            barcode_filtered.append(barcode)
        
    # 写出筛选的barcode文件，每个barcode写出一个文件
    # barcode_filtered_file = open(os.path(output_folder, f"{plot_name}_filter_barcode,.txt"),'w')
    # for k, v in barcode_count.items():
    #     if k in barcode_filtered:
    #         barcode_filtered_file.write("%s \t %s \n" % (k, v))
    # barcode_filtered_file.close()
    
    
    
    
    
    
    # 写出筛选的barcode对应的Sam文件
    

            
        
            
    
    
    
 
    
    

    
    
                
        
        
    
    
               
               
               
            
                
    
    
    


 