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
NK <- c('XCL1','KLRD1','KLRF1')
Tcells <- c('CD3D','CD3G','CD2')
Fibroblasts <- c('COL1A1','DCN','LUM')
Myeloids <- c('LYZ','CD68','TYROBP')
Epithelial <- c('CD24','KRT19','EPCAM')
Bcells <- c('CD79A','CD19','MS4A1')
Endothelial <- c('CLDN5','FLT1','RAMP2')
Plasma <- c('IGHG1','JCHAIN','MZB1')
DC <- c("LILRA4","CXCR3","IRF7")
Mast <- c("CPA3","TPSABT","TPSB2")

marker.list <- list(NK = NK,"T cell"=Tcells,'B cell'=Bcells,"Plasmas" =Plasma,
                    Myeloids=Myeloids,Fibroblasts=Fibroblasts,
                    Epithelials=Epithelial,Endothelials=Endothelial,
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
p
dev.off()

ggsave(filename='./FeaturePlot.png',height =20,width = 24,dpi=300, plot=p, bg='white')


# 绘制TOP10基因 热图
seurat.markers <- FindAllMarkers(scRNA_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
p_h = DoHeatmap(scRNA_seurat, features = top10$gene) + NoLegend() +
    theme(axis.text.y=element_text(size=8))
ggsave(filename = "/root/wangje/Project/吴霞/Data/01_Top10基因热图_new.png",plot=p_h, height=5,width=5,dpi=300)

# 写出表达矩阵
exp_mtx = as.matrix(scRNA_seurat@assays$RNA@counts) %>% as.data.frame()
exp_mtx$gene = rownames(exp_mtx)
data.table::fwrite(exp_mtx, file = "./exp_mtx.txt", sep="\t", row.names=T)

# 合并exp_matrix中的相同基因
df = pd.read_csv("./exp_mtx.txt",sep = '\t')
df_merge = df.groupby('gene'),sum()



gene_list = purrr::map(new_flist, function(x) rownames(x))
test = Reduce(union, gene_list)
expr_merge=aggregate(.~genes,sum,data=df)

####################################################################################
# 操作count
####################################################################################
suppressMessages(library(Matrix))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

combine_exon_intron <- function (df_gene, gene_count) 
{
    gene_count_exon = gene_count[df_gene$exon_intron == "exon", 
        ]
    gene_count_intron = gene_count[df_gene$exon_intron == "intron", 
        ]
    if (nrow(gene_count_exon) == nrow(gene_count_intron)) {
        gene_count_combine = gene_count_exon + gene_count_intron
    }
    else {
        gene_count_combine = gene_count_exon[-nrow(gene_count_exon), 
            ] + gene_count_intron
        gene_count_combine = rbind(gene_count_combine, gene_count_exon[nrow(gene_count_exon), 
            ])
    }
    return(gene_count_combine)
}

sciRNAseq_gene_count_summary <- function (gene_count_folder) {
    gene_matrix = paste(gene_count_folder, "/count.MM", sep = "")
    df_gene = paste(gene_count_folder, "/gene_name_annotate.txt", 
        sep = "")
    df_cell = paste(gene_count_folder, "/cell_annotate.txt", 
        sep = "")
    df_report = paste(gene_count_folder, "/report.MM", sep = "")
    report_annotate = paste(gene_count_folder, "/report_annotate.txt", 
        sep = "")
    df_gene = read.csv(df_gene, header = F)
    df_cell = read.csv(df_cell, header = F)
    gene_matrix = read.csv(gene_matrix, header = F)
    colnames(df_gene) = c("gene_id", "gene_type", "exon_intron", 
        "gene_name", "index")
    colnames(df_cell) = c("sample", "index")
    rownames(df_gene) = df_gene$gene_id
    rownames(df_cell) = df_cell$sample
    gene_count = sparseMatrix(i = gene_matrix$V1, j = gene_matrix$V2, 
        x = gene_matrix$V3)
    df_gene = df_gene[1:nrow(gene_count), ]
    rownames(gene_count) = df_gene$gene_id
    colnames(gene_count) = df_cell$sample
    gene_count = combine_exon_intron(df_gene, gene_count)
    df_gene = df_gene %>% filter(exon_intron == "exon")
    reportMM = read.csv(df_report, header = F)
    df_report = sparseMatrix(i = reportMM$V1, j = reportMM$V2, 
        x = reportMM$V3)
    df_report = as.matrix(t(df_report))
    df_report_annotate = read.csv(report_annotate, header = F)
    colnames(df_report) = df_report_annotate$V2
    df_report = data.frame(df_report)
    df_report["index"] = as.numeric(rownames(df_report))
    df_cell_combine = inner_join(df_cell, df_report, by = "index")
    df_cell_combine["all_exon"] = df_cell_combine$X.Perfect.intersect.exon.match + 
        df_cell_combine$X.Nearest.intersect.exon.match + df_cell_combine$X.Perfect.combine.exon.match + 
        df_cell_combine$X.Nearest.combine.exon.match
    df_cell_combine["all_intron"] = df_cell_combine$X.Perfect.intersect.gene.match + 
        df_cell_combine$X.Nearest.intersect.gene.match + df_cell_combine$X.Perfect.combine.gene.match + 
        df_cell_combine$X.Nearest.combine.gene.match
    df_cell_combine["all_reads"] = df_cell_combine$all_exon + 
        df_cell_combine$all_intron + df_cell_combine$X.No.match
    df_cell_combine["unmatched_rate"] = df_cell_combine$X.No.match/df_cell_combine$all_reads
    df_cell = df_cell_combine %>% select(sample, unmatched_rate)
    df_cell$UMI_count = df_cell_combine$all_exon + df_cell_combine$all_intron
    df_gene = df_gene %>% select(gene_id, gene_type, gene_name)
    return(list(df_cell, df_gene, gene_count))
}
result = sciRNAseq_gene_count_summary("./")
df_cell = result[[1]]
df_gene = result[[2]]
gene_count = result[[3]]
save(df_cell, df_gene, gene_count, file = paste0(output_folder, "/sci_summary.RData"))

###########################################################################################
library(data.table)
library(Seurat)
library(patchwork)
library(dplyr)


flist = list()
for(i in 1:96){
    input_file = paste0('/root/wangje/Project/吴霞/Data/07_count/BC230502-',i,'/','result4.txt')
    if(file.exists(input_file)){
        print(input_file)
        counts = fread(input_file, data.table = F,nThread = 10)
        rownames(counts) <- counts[, 1]
        counts <- counts[,-1]
        counts <- as(as.matrix(counts), "dgCMatrix")
        scRNA <- CreateSeuratObject(counts, simplify = T)
        scRNA[["percent.mt"]] <-PercentageFeatureSet(scRNA, pattern = "^MT-")
        flist[[i]] = scRNA
    }
}

# 合并flist文件
new_flist = Filter(Negate(is.null), flist)
flist = NULL 
scRNA = merge(new_flist[[1]], new_flist[2:length(new_flist)])

p1 = DimPlot(scRNA_seurat, group.by = 'RNA_snn_res.0.5') + 
    labs(title = '分辨率:RNA_snn_res.0.5',
        subtitle = paste0('genes: ',dim(scRNA_seurat)[1],'\t','cells: ',dim(scRNA_seurat)[2] ,'\n',
            '过滤条件:nFeature_RNA >= 50 & percent.mt < 80')
        )+
        theme(legend.position='bottom')

p2 = DimPlot(scRNA_seurat, group.by = 'celltype') + labs(title='SingleR注释')+theme(legend.position='bottom')



Tcells <- c('CD3D','CD3G','CD2')
Bcells<- c('CD79A','CD19','MS4A1')
NK <- c('XCL1','KLRD1','KLRF1')
Mast <- c('TPSAB1','TPSB2','CPA3','MS4A2')
DC <- c('CLEC4C','IL3RA','IRF7','LILRA4','CXCR3')
Fibroblasts <- c("COL1A1", "DCN", "LUM")
Myeloids <- c("LYZ", "CD68", "TYROBP")
Epithelial <- c("CD24", "KRT19", "EPCAM")
Endothelial <- c("CLDN5", "FLT1", "RAMP2")
Plasma <- c("IGHG1", "JCHAIN", "MZB1")

marker.list <- list("T cell" = Tcells, 
                    "B cell" = Bcells, 
                    NK = NK, "Plasmas" = Plasma, 
                    Myeloids = Myeloids,
                    Fibroblasts = Fibroblasts, 
                    Epithelials = Epithelial, 
                    Endothelials = Endothelial, 
                    Mast = Mast, 
                    DC = DC)


heatmap_color <- RColorBrewer::brewer.pal(name = "RdBu", n = 11)
pal <- rev(colorRampPalette(heatmap_color)(100))

Idents(scRNA_seurat) <- scRNA_seurat$celltype

p3 <- DotPlot(scRNA_seurat, assay = "RNA", dot.scale = 4, features = marker.list) +
    scale_color_gradientn(colours = pal) +
    annotate(geom = "segment", y = Inf, yend = Inf, color = "black", x = -Inf, xend = Inf, size = 1) +
    annotate(geom = "segment", x = Inf, xend = Inf, color = "black", y = -Inf, yend = Inf, size = 1) +
    annotate(geom = "segment", x = Inf, xend = Inf, color = "black", y = -Inf, yend = Inf, size = 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 12), strip.text.x = element_text(size = 10, angle = 360)) +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
    labs(x = "Gene Marker", y = "Cluster") + theme(legend.position = "right")

p5 = VlnPlot(scRNA_seurat, features = c('nFeature_RNA', 'nCount_RNA','percent.mt'), ncol = 3) & 
    theme(axis.title.x = element_blank())
design = 'ABCDEFFFFF
          ABCDEFFFFF'

p_all = p5 + p1+ p2 +p3  + plot_layout(design = design )

ggsave('./new_mt80.png', height = 5,width = 25, plot = p_all , bg = 'white')

ggsave('./Featureplot.png',height=12, width = 24, plot = p4 , bg = 'white')


flist = list()
for(i in list.files('./')){
    print(i)
    counts = fread(i, data.table = F,nThread = 10)
    rownames(counts) <- counts[, 1]
    counts <- counts[,-1]
    counts <- as(as.matrix(counts), "dgCMatrix")
    scRNA <- CreateSeuratObject(counts, simplify = T)
    scRNA[["percent.mt"]] <-PercentageFeatureSet(scRNA, pattern = "^MT-")
    flist[[i]] = scRNA

}


   ClusterID         celltype
1          0 Epithelial cells
2          1 Epithelial cells
3          2    T cells, CD8+
4          3 Epithelial cells
5          4 Epithelial cells
6          5 Epithelial cells
7          6      Fibroblasts
8          7    T cells, CD8+
9          8 Epithelial cells
10         9 Epithelial cells
11        10    Smooth muscle
12        11       Macrophage
13        12 Epithelial cells
14        13  Mesangial cells
15        14 Epithelial cells
16        15 Epithelial cells
17        16    T cells, CD4+
18        17  Mesangial cells
19        18  Mesangial cells
20        19 Epithelial_cells
21        20          B cells
22        21      Hepatocytes
23        22  Dendritic cells


scRNA_seurat@meta.data$celltype =  case_when(
    Idents(scRNA_seurat) %in% c(7,2,16) ~ 'T cells',
    Idents(scRNA_seurat) %in% c(6) ~ 'Fibroblasts',
    Idents(scRNA_seurat) %in% c(10) ~ 'Smooth muscle',
    Idents(scRNA_seurat) %in% c(20) ~ 'B cells',
    Idents(scRNA_seurat) %in% c(21) ~ 'Hepatocytes',
    Idents(scRNA_seurat) %in% c(22) ~ 'Dendritic cells',
    Idents(scRNA_seurat) %in% c(11) ~ 'Macrophage',
    TRUE ~ 'Epithelial cells'
)

genes = c('Col6a6','Glis1','Nr1h5','Col9a1','Ntng1','Trp63','Pth2r','Fndc3c1','Mybl1','Tfap2d',
          'C130060K24Rik','Dmbx1','Mylk4','Foxb1','Npy','Rab5a','Il31ra','Pax2','Id4','Emcn','Lamc3',
          'Tspan8','Mpz','Ppp1r1c','Cpa2','Hbb-bh1','Dlx6','Eomes','A1cf','Metrnl','Ms4a4a','Gmnc',
          'Uts2b','Myh6','Gp1ba','Tyr','Cryba2','Lcn2')
