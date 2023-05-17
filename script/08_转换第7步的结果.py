import sys
import os
import pandas as pd

"""首先将第7步生成的cell_annotate.txt，gene_name_annotate.txt和每个barcode的count文件转换为Seurat
   需要的barcodes.tsv，genes.tsv和matrix.mtx。"""

# 整理barcode文件   
def tidy_file(input_folder):
    if not os.path.exists(os.path.join(input_folder,"outfile")):
        os.mkdir(os.path.join(input_folder,"outfile"))
    else:
        print(f"输出文件已经存在: %s" %(os.path.join(input_folder,"outfile")))
    # 输入文件
    cell_input = os.path.join(input_folder,"cell_annotate.txt")
    gene_input = os.path.join(input_folder,"gene_name_annotate.txt")
    count_paths = [os.path.join(input_folder,file) for file in os.listdir(input_folder) if file.endswith(".count")]
    
    
    # 输出文件    
    cell_file = os.path.join(input_folder,"outfile","barcodes.tsv")
    gene_file = os.path.join(input_folder, "outfile","genes.tsv")
    exp_matrix = os.path.join(input_folder, "outfile","matrix_1.mtx")
    
    cell_count = []
    gene_count = []
    # 细胞文件
    print(f"{input_folder}:开始输出细胞文件！")
    with open(cell_input) as f_cell,open(cell_file, 'w') as fw:
        for line in f_cell.readlines():
            barcode = line.strip().split(",")[0]
            cell_count.append(barcode)
            fw.write(barcode + "\n")
    print("細胞文件结束")
    
    # 基因文件
    print(f"{input_folder}:开始输出基因文件！")
    with open(gene_input) as f_gene, open(gene_file , 'w') as fw:
        for line in f_gene.readlines():
            line = line.strip().split(',')
            gene_name = line[0]
            gene_symbol = line[-2]
            gene_count.append(gene_name)
            fw.write("%s \t %s \n" %(gene_name, gene_symbol))
    print("基因文件结束")
    
    # 表达文件
    print(f"{input_folder}:开始输出表达文件！")
    expression_pre = 0
    with open(exp_matrix,'w') as fw:
        for count_path in count_paths:
            with open(count_path, 'r') as fr:
                for line in fr.readlines():
                    expression_pre += 1
                    new_line = line.strip().replace(",", " ")
                    fw.write(new_line + "\n")
    df = pd.read_csv(os.path.join(input_folder,"outfile","matrix_1.mtx"),header=None)
    with open(os.path.join(input_folder,"outfile","matrix_1.mtx"), 'r') as fr, open(os.path.join(input_folder,"outfile","matrix.mtx"), 'w') as fw:
        fw.write("%%MatrixMarket matrix coordinate real general" + "\n")
        fw.write("%" + "\n")
        fw.write("%d %d %d\n" %(len(gene_count),len(cell_count),df.shape[0]))
        for line in fr.readlines():
            fw.write(line)
    os.remove(os.path.join(input_folder,"outfile","matrix_1.mtx"))
    df = None   
    
        

if __name__ == '__main__':
    input_folder = sys.argv[1]
    # input_folder = "/root/wangje/Project/吴霞/Data/06_splitSAM/BC230502-1/"
    tidy_file(input_folder=input_folder)
    
                
    
    
    
            
    