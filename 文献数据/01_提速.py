import os 
import sys
import gzip 
import time
import mappy as mp
import pandas as pd 

Read1 = "/root/wangje/Project/吴霞/SRR7827964/SRR7827964_1.fastq"
Read2 = "/root/wangje/Project/吴霞/SRR7827964/SRR7827964_2.fastq"
f3 = gzip.open("./outfile.fastq.gz", 'wb')

# ligation
ligiation = pd.read_csv("/root/wangje/Project/吴霞/SRR7827964/ligation_file.txt", header=None, names = ["ligation"])
ligiation_list = ligiation['ligation'].tolist()

# barcode
RT = pd.read_csv("/root/wangje/Project/吴霞/SRR7827964/RT_file.txt", header=None, names = ["RT"])
RT_list = RT['RT'].tolist()

match_read = {}
for name, seq, qual in mp.fastx_read(Read1):
    tmp_lig = seq[0:10]
    if tmp_lig in ligiation_list:
        tmp_RT = seq[len(tmp_lig) + 14 : len(tmp_lig) + 24]
        if tmp_RT in RT_list:
            umi = seq[len(tmp_lig)+6: len(tmp_lig) + 14]
            match_read[name] = "@" + name + "," + tmp_lig + tmp_RT + "," + umi

count = 0            
for name, seq, qual in mp.fastx_read(Read2):
    if name in match_read:
        count += 1
        # print(match_read[name])
        # print(seq)
        # print('+')
        # print(qual)
        f3.write((match_read[name] +"\n").encode())
        f3.write((seq + "\n").encode())
        f3.write(("+" + "\n").encode())
        f3.write((qual + "\n").encode())
        percent = (count/len(match_read))*100
        print('\r' + "#" * round(percent) + ">>> %d%% (%d/%d)" %(percent,count,len(match_read)),end="")
        
                
                        
            
            
            

