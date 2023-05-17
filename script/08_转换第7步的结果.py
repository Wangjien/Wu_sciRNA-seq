import sys
import os

"""首先将第7步生成的cell_annotate.txt，gene_name_annotate.txt和每个barcode的count文件转换为Seurat
   需要的barcodes.tsv，genes.tsv和matrix.mtx。"""

# 整理barcode文件   
def generate_barcode(input_folder, output_folder, sample_name):
    input_file = os.path.join(input_folder,sample_name)
    output_file = os.path.join(output_folder,"barcode.tsv")
    
    with open(input_file) as fr, open(output_file, 'w') as fw:
        for line in fr.readlines():
            # barcode = (line.strip().split(".")[1]).split(",")[0]
            barcode = line.strip().split(",")[0]
            fw.write(line)
    
# 整理genes文件            
def generate_genes(input_folder, output_folder, sample_name):
    input_file = os.path.join(input_folder,sample_name)
    output_file = os.path.join(output_folder, "genes.tsv")
    
    with open(input_file) as fr, open(output_file, 'w') as fw:
        for line in fr.readlines():
            line = line.strip().split(",")
            gene_name = line[0]
            gene_symbol = line[-2]
            fw.write("%s\t%s" %(gene_name, gene_symbol))

# 整理matrix文件
def generate_matrix(input_folder, output_folder):
    output_file = os.path.join(output_folder,"matrix.tsv")
    
    file_paths = [os.path.join(input_folder,file) for file in os.listdir(input_folder) if file.endswith(".count")]
    with open(output_file, 'w') as fw:
        for file_path in file_paths:
            print(file_path)
            with open(file_path) as fr:
                for line in fr.readlines():
                    new_line = line.rstrip().replace(",", " ")
                    fw.write(new_line + "\n")
                    # print(new_line)
                
                        
if __name__ == '__main__':
    
                
    
    
    
            
    