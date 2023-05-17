library(Seurat)
library(ggplot2)
library(dplyr)
library(purrr)
library(patchwork)
library(data.table)

# 生成Seurat文件,存储在列表中
flist = list()
for(i in seq_len(49)){
    pwd=paste0("/root/wangje/Project/吴霞/Data/06_splitSAM/","BC230502-",i,"/outfile")
    # print(pwd)
    if(dir.exists(pwd)){
        file_ok = file.exists(paste0(pwd,"/","barcodes.tsv")) & file.exists(paste0(pwd,"/","genes.tsv")) & file.exists(paste0(pwd,"/","matrix.mtx"))
        if(file_ok){
            data = Seurat::Read10X(data.dir = pwd)
            seurat_file = Seurat::CreateSeuratObject(counts = data, min.cells = 3) # min.features=200 去除了基因表达小于200的数据
            seurat_file[['percent,mt']] = PercentageFeatureSet(seurat_file, pattern = "^MT-")
            flist[[as.character(i)]] = seurat_file
        }else {
           cat(paste0("\033[31m", pwd,"缺少Seurat的輸入文件!", "\033[0m", "\n"))
           next
        }

    }
}
# 去除列表中的NULL
new_flist = Filter(Negate(is.null), flist)
flist = NULL 
# 绘制每个BC文件的小提琴图
plist = list()
for(i in seq_len(length(new_flist))){
    print(i)
    tryCatch({
    p = VlnPlot(new_flist[[i]], features=c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)&
        theme(axis.title.x = element_blank(),
              plot.title = element_text(size = 10))&
        labs(subtitle = paste0("BC230502-", names(new_flist[i])))
    plist[[i]] = p}, error = function(error) {
        cat(paste0("\033[31m",i,"有错误","\033[0m", "\n"))
        return(NULL)
    },finally={
        print("绘制结束！")

    })
}

# 保存图片
plist = Filter(Negate(is.null), plist) # 去除NULL
ggsave(filename = "/root/wangje/Project/吴霞/Data/小提琴图.png",plot=wrap_plots(plist,ncol=8),height=16,width=32,limitsize=FALSE)

# 合并列表中的Seurat对象为一个整体的Seurat对象
## 合并之前需要先去除dim为0的数据
new_flist <- lapply(new_flist, function(df) ifelse(dim(df)[1] > 0, df, NULL))
new_flist <- Filter(function(df) !is.null(df), new_flist)
scRNA_seurat = merge(new_flist[[1]],new_flist[2:length(new_flist)])

# 标准化及高变基因
scRNA_seurat <- NormalizeData(scRNA_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA_seurat <- FindVariableFeatures(scRNA_seurat, selection.method = "vst", nfeatures = 2000)

# 均一化
all.genes <- rownames(scRNA_seurat)
scRNA_seurat <- ScaleData(scRNA_seurat, features = all.genes)
# scRNA_seurat <- ScaleData(scRNA_seurat) # 默认是对高变基因进行归一化（默认2000个高变基因）

# pca
scRNA_seurat <- RunPCA(scRNA_seurat, features = VariableFeatures(object = scRNA_seurat))

# 聚类分群
scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:30)
scRNA_seurat <- FindClusters(scRNA_seurat, resolution = seq(0.05,1,0.05))
scRNA_seurat <- RunUMAP(scRNA_seurat, dims = 1:10)