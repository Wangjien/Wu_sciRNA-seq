suppressMessages(library(Matrix))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
args=commandArgs(T)

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

# 传入参数
result = sciRNAseq_gene_count_summary(args[1])
output_folder = args[2]
write.table(result[[1]],file = paste0(output_folder,result1.txt), sep = '\t', row.names=F, quote=F)


for(i in 24:89){
    pwd = paste0("/root/wangje/Project/吴霞/Data/06_splitSAM/",'BC230502-',i,'/outfile/')
    print(pwd)
    result = sciRNAseq_gene_count_summary(pwd)
    print(dim(result[[1]]))
    print(dim(result[[2]]))
    print(dim(result[[3]]))
    tryCatch({
    output_folder = paste0("/root/wangje/Project/吴霞/Data/07_count/","BC230502-",i)
    write.table(result[[1]], file = paste0(output_folder,'/','result1.txt'), sep = '\t', row.names=F, quote=F)
    write.table(result[[2]], file = paste0(output_folder,'/','result2.txt'), sep = '\t', row.names=F, quote=F)
    exp = as.data.frame(as.matrix(result[[3]]))
    exp$ID = rownames(exp)
    genefile = as.data.frame(result[[2]])
    exp$gene = genefile[match(exp$ID, genefile$gene_id), 'gene_name']
    exp = exp %>% dplyr::select('gene',everything())
    exp = exp[,-ncol(exp)]
    print(head(exp)[1:5,1:5])
    fwrite(exp, file = paste0(output_folder,'/','result3.txt'), sep = '\t', row.names=F, quote=F)
    },error=function(err){
        print('发生错误')
        print(err$message)
        next
    }, finally = {
        print('最终执行')
    }
    )
}


