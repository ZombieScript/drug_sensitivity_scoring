#################################### Step 0: 设置路径与加载包 ####################################
rm(list = ls())

# 确保已安装 apeglm
if (!requireNamespace("apeglm", quietly = TRUE)) {
  BiocManager::install("apeglm")
}

work_path <- "E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029"
setwd(work_path)

results_path <- file.path(work_path, "results")
rdata_path <- file.path(work_path, "Rdata")
dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
dir.create(rdata_path, showWarnings = FALSE, recursive = TRUE)

# 加载所需 R 包
library(DESeq2)
library(openxlsx)
library(org.Hs.eg.db)

#################################### Step 1: 加载表达矩阵与 metadata ####################################
# 读取表达矩阵
counts <- read.csv("counts.csv", header = TRUE, row.names = 1, check.names = FALSE)

# 读取 metadata 文件（此处指定编码为 GBK，如文件编码为 UTF-8，则改为 "UTF-8"）
metadata <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "GBK")

# 查看 metadata 前几行，确认文件结构
print(head(metadata))

# 清洗 SampleName 列（去除前后空格），确保与 counts 的列名一致
metadata$SampleName <- trimws(metadata$SampleName)

# 根据 Group 列筛选所需样本：只保留 parental 与 Resistant 两组样本
metadata_sub <- metadata[metadata$Group %in% c("parental", "Resistant"), ]

# 找出 counts 与 metadata_sub 中共有的样本
common_samples <- intersect(colnames(counts), metadata_sub$SampleName)
if (length(common_samples) == 0) {
  stop("未在 counts 文件与 metadata 文件中找到匹配的样本名，请检查文件内容！")
}

# 筛选交集样本
metadata_sub <- metadata_sub[metadata_sub$SampleName %in% common_samples, ]
counts_sub <- counts[, common_samples, drop = FALSE]

# 调整 metadata_sub 的顺序以匹配 counts_sub 的列顺序
metadata_sub <- metadata_sub[match(common_samples, metadata_sub$SampleName), ]
rownames(metadata_sub) <- metadata_sub$SampleName

# 将 counts_sub 数值取整（如有小数则四舍五入）
counts_sub <- round(counts_sub)

#################################### Step 2: 构建 DESeq2 数据对象 ####################################
dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                              colData = metadata_sub,
                              design = ~ Group)

# 设置对照组为 parental
dds$Group <- relevel(dds$Group, ref = "parental")

# 过滤低表达基因（总 reads 数至少为 10）
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 标准化并进行差异表达分析
dds <- DESeq(dds)

#################################### Step 3: 获取差异表达结果 ####################################
# 对比 Resistant vs parental，即 log2(Resistant/parental)
res <- results(dds, contrast = c("Group", "Resistant", "parental"))

# 使用 apeglm 方法收缩 logFC，提升结果稳定性
res <- lfcShrink(dds, coef = "Group_Resistant_vs_parental", res = res, type = "apeglm")

# 整理结果表，并按 padj 排序
res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj), ]

# 添加注释（基因全名）
res_df$GeneName <- mapIds(org.Hs.eg.db,
                          keys = rownames(res_df),
                          column = "GENENAME",
                          keytype = "SYMBOL",
                          multiVals = "first")

# 筛选显著差异基因（padj < 0.05 且 |log2FoldChange| >= 1）
filtered_res <- res_df[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1), ]
filtered_res <- filtered_res[order(filtered_res$log2FoldChange, decreasing = TRUE), ]

#################################### Step 4: 保存结果 ####################################
write.xlsx(res_df, file = file.path(results_path, "GSE200029_parental_vs_Resistant_all.xlsx"), rowNames = TRUE)
write.xlsx(filtered_res, file = file.path(results_path, "GSE200029_parental_vs_Resistant_filtered.xlsx"), rowNames = TRUE)

save(dds, res_df, filtered_res, file = file.path(rdata_path, "GSE200029_parental_vs_Resistant_DESeq2.RData"))

cat("✅ DESeq2 分析完成。\n")
cat("🔹 对照组: parental\n🔹 比较组: Resistant\n")
cat("🔹 差异基因数量: ", nrow(filtered_res), "\n")
cat("🔹 结果保存至: ", results_path, "\n")

