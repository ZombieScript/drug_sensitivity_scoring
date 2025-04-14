library(writexl)
library(dplyr)

cat("整合 DEG 与 Uniprot 并输出唯一基因-靶点数据...\n")

# 设置工作目录为 GSE20029 数据所在位置
setwd("E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029")

# 基于 ppi_final 提取 Uniprot-基因映射（双向提取 protein1 和 protein2 的注释）
gse20029_targets <- bind_rows(
  ppi_final %>% dplyr::select(GeneSymbol = protein1_gene, Target = protein1_uniprot),
  ppi_final %>% dplyr::select(GeneSymbol = protein2_gene, Target = protein2_uniprot)
) %>% 
  distinct(GeneSymbol, Target)

# 从 gse20029_results 中提取网络拓扑指标信息
# 注意：这里使用 gse20029_results 中的 gene_name 作为 GeneSymbol
gse20029_gene_topology <- gse20029_results %>%
  dplyr::select(GeneSymbol = gene_name, protein, Degree, Betweenness, Eigenvector, MCC, EPC)

# 内连接：整合 Uniprot 映射与拓扑信息（按 GeneSymbol 匹配），并去除重复的基因-靶点对
gse20029_final <- inner_join(gse20029_targets, gse20029_gene_topology, by = "GeneSymbol") %>%
  distinct(GeneSymbol, Target, .keep_all = TRUE)

# 添加拓扑权重打分（使用 gse20029_weights 对应的权重，注意权重名称为 "Degree_Shifted"、"Betweenness_Shifted" 等）
gse20029_final <- gse20029_final %>%
  mutate(
    Score = Degree * gse20029_weights["Degree_Shifted"] +
      Betweenness * gse20029_weights["Betweenness_Shifted"] +
      Eigenvector * gse20029_weights["Eigenvector_Shifted"] +
      MCC * gse20029_weights["MCC_Shifted"] +
      EPC * gse20029_weights["EPC_Shifted"]
  ) %>%
  arrange(desc(Score))

# 输出路径设置
output_path <- file.path("E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029", "results", "PPI")
output_file <- file.path(output_path, "Cleaned_GSE20029_Unique_Target_With_Topology.csv")

# 保存结果为 CSV 文件
write.csv(gse20029_final, output_file, row.names = FALSE)

cat("✅ 最终整合结果已保存至：\n", output_file, "\n")
cat("共输出唯一基因-靶点对数量：", nrow(gse20029_final), "\n")
