library(Seurat)
library(ggplot2)
library(dplyr)
library(purrr)
library(patchwork)
library(data.table)

# 生成Seurat文件,存储在列表中
flist = list()
input_folder = "/root/wangje/Project/吴霞/Data/06_splitSAM/"
for(i in seq_len(96)){
    pwd=paste0(input_folder,"BC230502-",i,"/outfile")
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
ggsave(filename = "/root/wangje/Project/吴霞/Data/小提琴图_分散_cutoff3.png",plot=wrap_plots(plist,ncol=12),height=24,width=48,limitsize=FALSE)

# 合并列表中的Seurat对象为一个整体的Seurat对象
## 合并之前需要先去除dim为0的数据
# new_flist <- lapply(new_flist, function(df) ifelse(dim(df)[1] > 0, df, NULL))
# new_flist <- Filter(function(df) !is.null(df), new_flist)
scRNA_seurat = merge(new_flist[[1]],new_flist[2:length(new_flist)])

# 查看合并后的数据分布(percent.mt, nFeature_RNA, nCount_RNA)
p = VlnPlot(scRNA_seurat, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, pt.size= 0.01)& theme(axis.title.x=element_blank())
ggsave(filename = "/root/wangje/Project/吴霞/Data/小提琴图_合并_cutoff3.png",plot=p,height=4,width=8,limitsize=FALSE)

# 进行筛选
scRNA_seurat_raw = scRNA_seurat
scRNA_seurat = scRNA_seurat[,scRNA_seurat$nFeature_RNA >= 200 & scRNA_seurat$percent.mt < 25]

# 标准化及高变基因
scRNA_seurat <- NormalizeData(scRNA_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA_seurat <- FindVariableFeatures(scRNA_seurat, selection.method = "vst", nfeatures = 10000)

# 均一化
all.genes <- rownames(scRNA_seurat)
scRNA_seurat <- ScaleData(scRNA_seurat, features = all.genes)
# scRNA_seurat <- ScaleData(scRNA_seurat) # 默认是对高变基因进行归一化（默认2000个高变基因）

# pca
scRNA_seurat <- RunPCA(scRNA_seurat, features = VariableFeatures(object = scRNA_seurat))

# 聚类分群
scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:30)
scRNA_seurat <- FindClusters(scRNA_seurat, resolution = seq(0.1,0.5,0.1))
scRNA_seurat <- RunUMAP(scRNA_seurat, dims = 1:30)

# 保存图片
ggsave(filename = "/root/wangje/Project/吴霞/Data/umap_new_cutoff3_保留nFeature200.png",plot=DimPlot(scRNA_seurat), height=5,width=5,dpi=300)
# ggsave(filename = "/root/wangje/Project/吴霞/Data/Vlnplot_整体.png",plot=VlnPlot(scRNA_seurat, features=c('nCount_RNA','nFeature_RNA','percent.mt'), ncol=3), height=5,width=8,dpi=300)

# 绘制FeaturePlot
NKT <- c('CD3D','CD3G','CD2')
Fibroblasts <- c('COL1A1','DCN','LUM')
Myeloids <- c('LYZ','CD68','TYROBP')
Epithelial <- c('CD24','KRT19','EPCAM')
Bcells <- c('CD79A','CD19','MS4A1')
Endothelial <- c('CLDN5','FLT1','RAMP2')
Plasma <- c('IGHG1','JCHAIN','MZB1')
Hepatocytes <- c('ALB','APOB','HP')
Keratinocytes <- c("KRT5","KRT14","FABP5")
DC <- c("LILRA4","CXCR3","IRF7")
Mast <- c("CPA3","TPSABT","TPSB2")

marker.list <- list("NK&T cell"=NKT,'B cell'=Bcells,"Plasmas" =Plasma,
                    Myeloids=Myeloids,Fibroblasts=Fibroblasts,
                    Epithelials=Epithelial,Endothelials=Endothelial,
                    Hepatocytes=Hepatocytes,Keratinocytes=Keratinocytes,

                    DC = DC, Mast = Mast)

plotFeature <- function(scRNA_data = scRNA_data,
                        choose = "Feature",
                        col_num = 6, marker.list = marker.list,...) {
    pacman::p_load("Seurat", "ggplot2", "tidyverse")
    DefaultAssay(scRNA_data) <- "RNA"
    plist <- list()
    if (is.null(choose)) {
        message("请选择绘图类型")
    } else if (choose == "Feature") {
        for (i in names(marker.list)) {
            for (j in marker.list[[i]]) {
                #    print(paste0(i,"_",j))
                tmp <- tryCatch(
                    {
                        FeaturePlot(scRNA_data, features = j) +
                            theme(legend.position = "right") +
                            labs(title = paste0(i, "_", j))
                    },
                    error = function(e) {
                        message("Error @ ", j)
                        return(NA)
                    },
                    finally = {
                        message(paste0(i, "_", j, "_next..."))
                    }
                )
                plist[[paste0(i, "_", j)]] <- tmp
            }
        }
        p_new <- Filter(Negate(anyNA), plist)
        p <- wrap_plots(p_new, bycol = T, ncol = col_num)
        return(p)}}
png("/root/wangje/Project/吴霞/Data/大群markerFeaturePlot_cutoff3_new.png",height =2000,width = 5000,res=300)
png("/root/wangje/Project/吴霞/Data/GSM4186980.png",height =2000,width = 5000,res=300)
p = plotFeature(scRNA_data=test,choose="Feature",col_num=6,marker.list=marker.list)
dev.off()

ggsave(filename='./FeaturePlot.png',height =20,width = 24,dpi=300, plot=p, bg='white')


# 绘制TOP10基因 热图
seurat.markers <- FindAllMarkers(scRNA_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
p_h = DoHeatmap(scRNA_seurat, features = top10$gene) + NoLegend() +
    theme(axis.text.y=element_text(size=3))
ggsave(filename = "/root/wangje/Project/吴霞/Data/01_Top10基因热图_new.png",plot=p_h, height=5,width=5,dpi=300)

# 写出表达矩阵
exp_mtx = as.matrix(scRNA_seurat@assays$RNA@counts) %>% as.data.frame()
exp_mtx$gene = rownames(exp_mtx)
data.table::fwrite(exp_mtx, file = "./exp_mtx.txt", sep="\t", row.names=T)

# 合并exp_matrix中的相同基因
df = pd.read_csv("./exp_mtx.txt",sep = '\t')



gene_list = purrr::map(new_flist, function(x) rownames(x))
test = Reduce(union, gene_list)
expr_merge=aggregate(.~genes,sum,data=df)

