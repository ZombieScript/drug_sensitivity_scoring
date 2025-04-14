###############################################################################
# 1. 清除环境变量
###############################################################################
rm(list = ls())

###############################################################################
# 2. 加载必要库
###############################################################################
library(dplyr)
library(tidyr)
library(readxl)

###############################################################################
# 3. 读取主要数据文件路径 (根据实际路径进行修改)
###############################################################################
gse20029_file  <- "E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029/results/PPI/Cleaned_GSE20029_Unique_Target_With_Topology.csv"
drug_file      <- "E:/Bioinformatics/download/ChEMBL/DrugDbPathway/陶术抗肿瘤药列表/drug_target_mechanism_final.csv"
topo_file      <- "E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029/results/PPI/TopoNetwork_GSE20029.RData"
deg_excel_file <- "E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029/results/GSE200029_parental_vs_Resistant_filtered.xlsx"

###############################################################################
# 4. 加载数据
###############################################################################
cat("加载去重后的 PPI 数据...\n")
gse20029_data <- read.csv(gse20029_file, header = TRUE, stringsAsFactors = FALSE)

cat("加载药物数据库信息...\n")
drug_data  <- read.csv(drug_file, header = TRUE, stringsAsFactors = FALSE)

cat("加载拓扑网络权重数据...\n")
load(topo_file)  # 假设里面有 gse20029_results, gse20029_weights, weight_df 等

###############################################################################
# 5. 从 Excel 中读取 log2FoldChange 信息，准备合并到 gse20029_results
###############################################################################
cat("从 Excel 中读取 GSE200029 差异表达 (log2FoldChange) 信息...\n")
deg_data <- read_excel(deg_excel_file, sheet = 1)
colnames(deg_data)[1] <- "gene_name"
cat("Excel 列名:\n")
print(colnames(deg_data))

###############################################################################
# 6. 将 log2FoldChange 合并到 gse20029_results
###############################################################################
cat("将 log2FoldChange 合并进 gse20029_results...\n")
gse20029_results <- gse20029_results %>%
  left_join(deg_data[, c("gene_name", "log2FoldChange")], by = "gene_name")

cat("\ngse20029_results 合并后列名:\n")
print(colnames(gse20029_results))

###############################################################################
# 7. 确认 gse20029_results 中是否含有必须列: 'gene_name', 'uniprot_id'
###############################################################################
required_cols <- c("gene_name", "uniprot_id")
if(!all(required_cols %in% colnames(gse20029_results))){
  stop("gse20029_results 数据中缺少 'gene_name' 或 'uniprot_id' 列。请确认 RData 内容。")
}

###############################################################################
# 8. 将 gse20029_data 与 gse20029_results 合并以获取 uniprot_id、Shifted 指标 以及 log2FoldChange
###############################################################################
cat("\n合并 gse20029_data 与 gse20029_results...\n")
gse20029_merged <- merge(
  gse20029_data,
  gse20029_results[, c("gene_name", "uniprot_id",
                       "Degree_Shifted", "Betweenness_Shifted", "Eigenvector_Shifted",
                       "MCC_Shifted", "EPC_Shifted",
                       "log2FoldChange")],
  by.x = "GeneSymbol",
  by.y = "gene_name",
  all.x = TRUE
)

###############################################################################
# 9. 检查合并后的 gse20029_merged 是否包含 'uniprot_id'
###############################################################################
cat("\n检查合并后的数据是否包含 'uniprot_id'...\n")
if(!"uniprot_id" %in% colnames(gse20029_merged)){
  stop("合并后数据中缺少 'uniprot_id' 列。请检查合并步骤。")
}

###############################################################################
# 10. 清洗 Uniprot ID 格式 (大写、去空格、截取前6位等)
###############################################################################
cat("\n清洗 Uniprot ID 的格式...\n")
gse20029_merged <- gse20029_merged %>%
  mutate(
    uniprot_id = toupper(trimws(uniprot_id)),
    uniprot_id = substr(uniprot_id, 1, 6)
  )

drug_data <- drug_data %>%
  mutate(
    uniprot_id = toupper(trimws(uniprot_id)),
    uniprot_id = substr(uniprot_id, 1, 6)
  )

###############################################################################
# 11. 检查 Uniprot ID 重叠情况
###############################################################################
cat("\n检查 Uniprot ID 的重叠情况...\n")
unique_gse20029_uniprot <- unique(gse20029_merged$uniprot_id)
unique_drug_uniprot     <- unique(drug_data$uniprot_id)
overlap_count           <- length(intersect(unique_gse20029_uniprot, unique_drug_uniprot))
cat("GSE20029 与 Drug 数据集的 Uniprot ID 重叠数量: ", overlap_count, "\n")

###############################################################################
# 12. 计算 Target_importance (加权总和)
###############################################################################
cat("\n根据 weight_df 为各拓扑指标赋权...\n")
GSE20029_weight_df <- weight_df
method_map <- c(
  "Degree_Zscore"      = "Degree_Z",
  "Betweenness_Zscore" = "Betweenness_Z",
  "Eigenvector_Zscore" = "Eigenvector_Z",
  "MCC_Zscore"         = "MCC_Z",
  "EPC_Zscore"         = "EPC_Z"
)
norm_method_names <- method_map[GSE20029_weight_df$Method]
method_weights <- setNames(GSE20029_weight_df$GSE20029_Weight, norm_method_names)

cat("\n对 gse20029_merged 计算 Target_importance...\n")
gse20029_merged <- gse20029_merged %>%
  rowwise() %>%
  mutate(
    Target_importance =
      (Degree_Shifted      * method_weights["Degree_Z"]) +
      (Betweenness_Shifted * method_weights["Betweenness_Z"]) +
      (Eigenvector_Shifted * method_weights["Eigenvector_Z"]) +
      (MCC_Shifted         * method_weights["MCC_Z"]) +
      (EPC_Shifted         * method_weights["EPC_Z"])
  ) %>%
  ungroup()

###############################################################################
# 13. (可选) 生成 Adjusted_logFC
###############################################################################
if("log2FoldChange" %in% colnames(gse20029_merged)){
  cat("\n生成 Adjusted_logFC...\n")
  gse20029_merged <- gse20029_merged %>%
    mutate(Adjusted_logFC = sqrt(abs(log2FoldChange)) * sign(log2FoldChange))
} else {
  cat("提示: gse20029_merged 中没有 'log2FoldChange'，无法生成 Adjusted_logFC。\n")
}

###############################################################################
# 14. (如上已处理 Adjusted_logFC，则此处可省略或保留)
###############################################################################

###############################################################################
# 15. 去重 (distinct)
###############################################################################
cat("\n去除重复行...\n")
gse20029_data_clean <- gse20029_merged %>%
  distinct(GeneSymbol, Target, .keep_all = TRUE)

cat("去重后 gse20029_data_clean 的行数: ", nrow(gse20029_data_clean), "\n")

###############################################################################
# 16. 再次清洗 Target 列 (如需要截取前6位)
###############################################################################
cat("\n清洗 Target 列...\n")
gse20029_data_clean <- gse20029_data_clean %>%
  mutate(
    Target = toupper(trimws(Target)),
    Target = substr(Target, 1, 6)
  )

###############################################################################
# 17. 合并药物信息 (left_join)
###############################################################################
cat("\n合并药物信息...\n")
gse20029_data_final <- gse20029_data_clean %>%
  left_join(drug_data, by = c("Target" = "uniprot_id"))

###############################################################################
# 18. 检查合并情况
###############################################################################
cat("\n检查合并后的数据...\n")
cat("GSE20029 数据中匹配到药物信息的数量: ", sum(!is.na(gse20029_data_final$drug_name)), "\n")

###############################################################################
# 19. 过滤掉未匹配到药物信息的记录
###############################################################################
cat("\n过滤掉未匹配的药物信息...\n")
gse20029_matched <- gse20029_data_final %>%
  filter(!is.na(drug_name))

###############################################################################
# 20. 计算积分 (Total_Score)
###############################################################################
calculate_total_score <- function(df, method_weights) {
  df %>%
    rowwise() %>%
    mutate(
      Degree_Score      = log2FoldChange * ((1 / standard_value)^(1/2)) * Degree_Shifted      * 100 * method_weights["Degree_Z"],
      Betweenness_Score = log2FoldChange * ((1 / standard_value)^(1/2)) * Betweenness_Shifted * 100 * method_weights["Betweenness_Z"],
      Eigenvector_Score = log2FoldChange * ((1 / standard_value)^(1/2)) * Eigenvector_Shifted * 100 * method_weights["Eigenvector_Z"],
      MCC_Score         = log2FoldChange * ((1 / standard_value)^(1/2)) * MCC_Shifted         * 100 * method_weights["MCC_Z"],
      EPC_Score         = log2FoldChange * ((1 / standard_value)^(1/2)) * EPC_Shifted         * 100 * method_weights["EPC_Z"],
      Raw_Total_Score = sum(
        c(Degree_Score, Betweenness_Score, Eigenvector_Score, MCC_Score, EPC_Score),
        na.rm = TRUE
      ),
      effect_multiplier = case_when(
        effect_type == "suppress" ~ 1,
        effect_type == "promote"  ~ -1,
        TRUE                      ~ 1
      ),
      Total_Score = Raw_Total_Score * effect_multiplier
    ) %>%
    ungroup() %>%
    dplyr::select(-Raw_Total_Score, -effect_multiplier)
}

cat("\n计算 gse20029_matched 的积分...\n")
gse20029_scored <- calculate_total_score(gse20029_matched, method_weights)

###############################################################################
# 21. (可选) 如果需要统一列名
###############################################################################
# colnames(gse20029_scored)[colnames(gse20029_scored) == "GeneSymbol"] <- "GSE20029_GeneSymbol"

###############################################################################
# 22. 按药物进行汇总打分
###############################################################################
cat("\n对每种药物进行汇总统计...\n")
drug_scores_gse20029 <- gse20029_scored %>%
  group_by(drug_name) %>%
  summarise(
    Total_Score           = sum(Total_Score, na.rm = TRUE),
    Total_Importance      = sum(Target_importance, na.rm = TRUE),
    Matched_Targets       = n(),
    Average_StandardValue = mean(standard_value, na.rm = TRUE),
    Matched_Gene_Names    = paste(unique(GeneSymbol), collapse = "; "),
    Matched_Target_Details = paste0(
      "Gene: ", GeneSymbol,
      ", Target: ", uniprot_id,
      ", Standard Value: ", standard_value,
      ", logFC: ", log2FoldChange,
      ", Effect Type: ", effect_type
    ) %>% paste(collapse = "; "),
    .groups = "drop"
  )

###############################################################################
# 23. 查看结果示例
###############################################################################
cat("\nGSE20029 药物汇总打分示例:\n")
print(head(drug_scores_gse20029))

###############################################################################
# 24. 保存结果
###############################################################################
cat("\n保存药物汇总打分...\n")
write.csv(drug_scores_gse20029,
          "E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029/results/drug_scores_gse20029.csv",
          row.names = FALSE)

save(gse20029_scored,
     file = "E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029/results/DrugMappingResults_GSE20029.RData")

###############################################################################
# 25. 输出靶点重要性文件
###############################################################################
output_dir <- "E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029/results/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("\n保存 GSE20029 靶点重要性...\n")
write.csv(
  gse20029_data_clean %>%
    dplyr::select(GeneSymbol, uniprot_id, log2FoldChange, Target_importance, everything()),
  file = paste0(output_dir, "gse20029_target_importance.csv"),
  row.names = FALSE
)

###############################################################################
# 26. (可选) 如需 RData 格式保存更多对象
###############################################################################
rdata_file <- "E:/Bioinformatics/Method_Validation/GSE200029/Rdata/"
if(!dir.exists(rdata_file)) dir.create(rdata_file, recursive = TRUE)
save(gse20029_scored, gse20029_data_clean,
     file = paste0(rdata_file, "TargetImportanceResults_GSE20029.RData"))

###############################################################################
# 27. 结束
###############################################################################
cat("\n处理完成。\n")

