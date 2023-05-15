import itertools
import sys
import argparse
import pandas as pd

def generate_combined_barcode(ligation_file, RT_barcode_file, p7_file):
    ligation = pd.read_csv(ligation_file, header=None, names=['ligation'], sep = '\t')
    ligation_list = ligation["ligation"].tolist()
    
    RT_barcode = pd.read_csv(RT_barcode_file, header=None, names=['RT'], sep ='\t')
    RT_list = RT_barcode["RT"].tolist()
    
    p7 = pd.read_csv(p7_file, header=None, names=['sampleID','p7'], sep = '\t')
    p7_barcode_list = p7["p7"].tolist()
    
    combined = [''.join(item) for item in itertools.product(RT_list, ligation_list, p7_barcode_list)]
    for item in combined:
        print(item)

if __name__ == "__main__":
    generate_combined_barcode(
        ligation_file= sys.argv[1],
        RT_barcode_file= sys.argv[2],
        p7_file= sys.argv[3]
    )
     
# 用法 python 生成barcode的所有组合文件.py ligation_barcode.txt RT_barcode.txt sample_p7.txt > combined_barcode.txt    
