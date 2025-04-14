#################################### Step 0: 清理环境与加载包 ####################################
rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(openxlsx)  # 读取 Excel 文件

# 设置工作目录（更新为 GSE200029 的目录）
setwd("E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029")

# 定义路径
alias_file <- "E:/Bioinformatics/download/STRING/9606.protein.aliases.v12.0.txt.gz"
links_file <- "E:/Bioinformatics/download/STRING/9606.protein.links.v12.0.txt.gz"
deg_file   <- "results/GSE200029_parental_vs_Resistant_filtered.xlsx"

output_path <- "results/PPI"
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

#################################### Step 1: 读取 DEG 基因列表 ####################################
message("读取差异表达基因列表 ...")
deg_data <- openxlsx::read.xlsx(deg_file, colNames = TRUE)

# 手动设置列名（根据文件列顺序）
colnames(deg_data) <- c("Gene", "log2FoldChange", "pvalue", "padj", "GeneName")

# 确保存在 Gene 列
if (!("Gene" %in% colnames(deg_data))) {
  stop("差异表达结果中未包含 Gene 列，请检查输入文件格式。")
}

deg_genes <- unique(deg_data$Gene)
cat("显著差异基因数（去重后）:", length(deg_genes), "\n")

#################################### Step 2: 读取 STRING alias 映射表 ####################################
message("读取 STRING alias 映射表 ...")
alias_data <- fread(alias_file, sep = "\t", header = TRUE, showProgress = TRUE)

message("整理 gene_symbol 和 UniProt ID 映射 ...")

# 提取 gene symbol 映射
alias_gene <- alias_data[source %in% c("Ensembl_HGNC_symbol", "Ensembl_HGNC"),
                         .(string_protein_id = `#string_protein_id`, gene_symbol = alias)]

alias_gene_unique <- alias_gene[, .SD[1], by = string_protein_id]

# 提取 UniProt 映射（格式限制）
alias_uniprot <- alias_data[grepl("^[OPQ][0-9A-Z]{4}[0-9]$", alias) & grepl("UniProt", source),
                            .(string_protein_id = `#string_protein_id`, uniprot_id = alias)]

alias_uniprot_unique <- alias_uniprot[, .SD[1], by = string_protein_id]

# 构建映射表
mapping_table <- merge(alias_gene_unique, alias_uniprot_unique, by = "string_protein_id", all.x = TRUE)

# 释放内存
rm(alias_data, alias_gene, alias_uniprot, alias_gene_unique, alias_uniprot_unique)
gc()

#################################### Step 3: 映射 DEG 到 STRING ID ####################################
message("映射差异基因至 STRING ID ...")
deg_proteins <- mapping_table[gene_symbol %in% deg_genes,
                              .(string_protein_id, gene_symbol, uniprot_id)]
deg_protein_ids <- unique(deg_proteins$string_protein_id)

cat("成功映射 STRING 蛋白数:", length(deg_protein_ids), "\n")

#################################### Step 4: 读取 STRING PPI 数据 ####################################
message("读取 PPI 数据 ...")
ppi_data <- fread(links_file, sep = " ", header = TRUE, showProgress = TRUE)

message("筛选高可信度互作 (score ≥ 700) 且在 DEG 蛋白集中 ...")
ppi_filtered <- ppi_data[combined_score >= 700 &
                           protein1 %in% deg_protein_ids &
                           protein2 %in% deg_protein_ids]

cat("筛选后 PPI 记录数:", nrow(ppi_filtered), "\n")

rm(ppi_data)
gc()

#################################### Step 5: 注释基因符号和 UniProt ####################################
message("注释 PPI 网络 ...")
setkey(mapping_table, string_protein_id)

# 注释 protein1
ppi_filtered <- merge(ppi_filtered,
                      mapping_table[, .(string_protein_id, gene_symbol, uniprot_id)],
                      by.x = "protein1",
                      by.y = "string_protein_id",
                      all.x = TRUE)
setnames(ppi_filtered, c("gene_symbol", "uniprot_id"), c("protein1_gene", "protein1_uniprot"))

# 注释 protein2
ppi_filtered <- merge(ppi_filtered,
                      mapping_table[, .(string_protein_id, gene_symbol, uniprot_id)],
                      by.x = "protein2",
                      by.y = "string_protein_id",
                      all.x = TRUE)
setnames(ppi_filtered, c("gene_symbol", "uniprot_id"), c("protein2_gene", "protein2_uniprot"))

#################################### Step 6: 去除方向性重复 ####################################
message("去除方向性重复 ...")
ppi_filtered[, pair := ifelse(protein1 < protein2,
                              paste(protein1, protein2, sep = "_"),
                              paste(protein2, protein1, sep = "_"))]

ppi_unique <- unique(ppi_filtered, by = "pair")

ppi_final <- ppi_unique[, .(protein1, protein2, combined_score,
                            protein1_gene, protein1_uniprot,
                            protein2_gene, protein2_uniprot)]

cat("去重后 PPI 网络大小:", nrow(ppi_final), "\n")

#################################### Step 7: 保存结果 ####################################
message("保存 PPI 网络结果 ...")

# 保存 TSV 文件
fwrite(ppi_final,
       file = file.path(output_path, "GSE200029_HighConfidence_PPI_With_Gene_and_UniProt.tsv"),
       sep = "\t")

# 保存 RData 文件
save(ppi_final, deg_proteins,
     file = file.path(output_path, "GSE200029_PPIAnalysis_highconfidence.RData"))

cat("✅ GSE200029 PPI 网络分析完成！\n",
    "结果保存路径: ", output_path, "\n")
