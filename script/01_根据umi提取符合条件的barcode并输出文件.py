import os
import sys
import gzip
import time
import pandas as pd
import pyfastx
from multiprocessing import pool

def Extract_barcode(input_folder, output_folder, sample_id, ligation_barcode_file, RT_Barcode_file, p7_file):
    """
    从Read1序列中根据umi序列的位置,提取RT barcode和ligation barcode,添加到对应的双端测序的Read2序列,并输出结果文件。
    """
    time.start = time.time()
    f1 = os.path.join(input_folder,f'{sample_id}_R1.fq.gz')
    f2 = os.path.join(input_folder,f'{sample_id}_R2.fq.gz')
    f3 = os.path.join(output_folder,f"{sample_id}_R2_barcode.fq.gz")
    f3_w = gzip.open(f3, 'wb')
    
    # 读入ligation barcode和RT barcode和p7的文件文件
    ligation_file = pd.read_csv(ligation_barcode_file, header=None, names=['ligation'],sep='\t')
    ligation_list = ligation_file['ligation'].tolist()
    RT_file = pd.read_csv(RT_Barcode_file, header=None, names=['RT_barcode'],sep='\t')
    RT_list = RT_file['RT_barcode'].to_list()
    f1 = os.path.join(input_folder,f'{sample_id}_R1.fq.gz')
    f2 = os.path.join(input_folder,f'{sample_id}_R2.fq.gz')
    f3 = os.path.join(output_folder,f"{sample_id}_R2_barcode.fq.gz")
    f3_w = gzip.open(f3, 'wb')
    log = open(os.path.join(output_folder,f"{sample_id}_log.txt"),'w')
    
    # 读入ligation barcode和RT barcode和p7的文件文件
    ligation_file = pd.read_csv(ligation_barcode_file, header=None, names=['ligation'],sep='\t')
    ligation_list = ligation_file['ligation'].tolist()
    RT_file = pd.read_csv(RT_Barcode_file, header=None, names=['RT_barcode'],sep='\t')
    RT_list = RT_file['RT_barcode'].to_list()
    p7 = pd.read_csv(p7_file, header=None, names=['sample','p7_barcode'],sep='\t')
    p7_dict = dict(zip(p7['sample'].tolist(),p7['p7_barcode'].tolist()))
    
    # 接受符合条件的序列
    ligation_match = []
    RT_match=[]
    name_match = []
    umi_match = []
    for name,seq,qual in pyfastx.Fastq(f1, build_index=False):
        if 'CAGAGC' in seq:
            umi_index = seq.find('CAGAGC')
            match_ligation = seq[int(umi_index -10):int(umi_index -1)]
            match_RT = seq[int(umi_index + 6+8):int(umi_index+6+18)]
            if match_RT in RT_list and match_ligation in ligation_list:
                ligation_match.append(match_ligation)
                RT_match.append(match_RT)
                umi = seq[int(umi_index):int(umi_index+14)]
                name_match.append(name)
                umi_match.append(umi)
    merge_barcode_umi = {k:v1+v2+str(p7_dict[sample_id])+','+umi_seq for k,v1,v2, umi_seq in zip(name_match,RT_match,ligation_match,umi_match)}
    # 从Read2序列中查找符合的条件
    count = 0
    reads_count = []
    for name,seq,qual in pyfastx.Fastq(f2,build_index=False):  
        reads_count.append(name)
        
        if name in merge_barcode_umi:
            # print('@' +name + ',' +  merge_barcode_umi[name])
            # print(seq)
            # print('+')
            # print(qual)
            f3_w.write(('@'+name + ',' + merge_barcode_umi[name]+'\n').encode())
            f3_w.write((seq + '\n').encode())
            f3_w.write(('+' + '\n').encode())
            f3_w.write((qual + '\n').encode())
            count +=1
            percent = (count/len(merge_barcode_umi))*100
            print('\r' + ">>>>>>>>>>>>>>>>>>>>>>> %d%% (%d/%d)" %(percent,count,len(merge_barcode_umi)),end="")
            # print(count)
    time.end = time.time()
    log.write(
        """
        文件:%s
        测序总reads数目:%s,
        符合条件的reads数目:%s,
        占比:%.4f,
        耗时:%s 秒
        """ % (sample_id,len(reads_count),len(merge_barcode_umi),
               len(merge_barcode_umi)/len(reads_count),
               time.end - time.start))
    
if __name__ == '__main__':
    # input_folder = '/home/data/vip9t22/test/BC230502-1'
    # output_folder = '/home/data/vip9t22/test/BC230502-1'
    # sample_id = 'BC230502-1'
    # ligation_barcode_file = '/home/data/vip9t22/test/BC230502-1/ligation_barcode.txt'
    # RT_Barcode_file = '/home/data/vip9t22/test/BC230502-1/RT_barcode.txt'
    # p7_file = '/home/data/vip9t22/test/BC230502-1/sample_p7.txt'
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    sample_id = sys.argv[3]
    ligation_barcode_file = sys.argv[4]
    RT_Barcode_file = sys.argv[5]
    p7_file = sys.argv[6]
    Extract_barcode(
        input_folder = input_folder,
        output_folder = output_folder,
        sample_id = sample_id,
        ligation_barcode_file = ligation_barcode_file,
        RT_Barcode_file = RT_Barcode_file,
        p7_file = p7_file
    )





