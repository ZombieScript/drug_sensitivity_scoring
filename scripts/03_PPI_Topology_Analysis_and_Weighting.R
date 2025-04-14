# ==============================================
# GSE20029 差异表达基因对应 PPI 网络的拓扑分析
# ==============================================

rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(igraph)
library(ggplot2)
library(reshape2)

# 设置工作目录（更新为 GSE200029 的目录）
work_path <- "E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/GSE200029"
setwd(work_path)

# 加载 PPI 网络数据（使用之前保存的 GSE20029 对象）
ppi_data_path <- file.path(work_path, "results", "PPI", "GSE200029_PPIAnalysis_highconfidence.RData")
if (!file.exists(ppi_data_path)) {
  stop("错误：未找到 PPI 网络数据文件: ", ppi_data_path)
}
load(ppi_data_path)  # 加载对象： ppi_final, deg_proteins

# 使用加载的 ppi_final 对象
ppi_unique <- ppi_final

# 构建 igraph 网络
gse20029_network <- graph_from_data_frame(ppi_unique[, c("protein1", "protein2")], directed = FALSE)

# 节点注释函数
create_node_attributes <- function(ppi_data) {
  node1 <- unique(ppi_data[, c("protein1", "protein1_gene", "protein1_uniprot")])
  colnames(node1) <- c("protein", "gene_name", "uniprot_id")
  node2 <- unique(ppi_data[, c("protein2", "protein2_gene", "protein2_uniprot")])
  colnames(node2) <- c("protein", "gene_name", "uniprot_id")
  
  node_all <- bind_rows(node1, node2) %>%
    group_by(protein) %>%
    summarise(
      gene_name = paste(unique(gene_name), collapse = "; "),
      uniprot_id = paste(unique(uniprot_id), collapse = "; "),
      .groups = "drop"
    )
  return(node_all)
}

gse20029_node_attributes <- create_node_attributes(ppi_unique)

# Z-score 平移函数（使所有分值变为正值）
shift_to_positive <- function(z_scores, epsilon = 0.01) {
  min_z <- min(z_scores, na.rm = TRUE)
  c <- epsilon - min_z
  return(z_scores + c)
}

# 拓扑指标计算函数
compute_topology_metrics <- function(network, node_attributes) {
  results <- data.frame(protein = V(network)$name, stringsAsFactors = FALSE)
  
  results$Degree <- degree(network)
  results$Betweenness <- betweenness(network)
  results$Eigenvector <- eigen_centrality(network)$vector
  
  results$MCC <- sapply(V(network), function(node) {
    neighbors_node <- ego(network, order = 1, nodes = node)[[1]]
    subgraph <- induced_subgraph(network, neighbors_node)
    cliques_sub <- cliques(subgraph)
    if (length(cliques_sub) == 0) return(0)
    return(max(sapply(cliques_sub, length)))
  })
  
  edge_bet <- edge_betweenness(network)
  results$EPC <- sapply(V(network), function(node) {
    incident_edge_indices <- incident(network, node, mode = "all")
    sum(edge_bet[incident_edge_indices], na.rm = TRUE)
  })
  
  # 注释节点的基因信息
  results <- left_join(results, node_attributes, by = "protein")
  
  # 计算 Z-score，并平移至正数
  results <- results %>%
    mutate(
      Degree_Z = scale(Degree)[,1],
      Betweenness_Z = scale(Betweenness)[,1],
      Eigenvector_Z = scale(Eigenvector)[,1],
      MCC_Z = scale(MCC)[,1],
      EPC_Z = scale(EPC)[,1],
      
      Degree_Shifted = shift_to_positive(Degree_Z),
      Betweenness_Shifted = shift_to_positive(Betweenness_Z),
      Eigenvector_Shifted = shift_to_positive(Eigenvector_Z),
      MCC_Shifted = shift_to_positive(MCC_Z),
      EPC_Shifted = shift_to_positive(EPC_Z)
    )
  
  return(results)
}

# 计算网络的拓扑指标
gse20029_results <- compute_topology_metrics(gse20029_network, gse20029_node_attributes)
head(gse20029_results)

# 检查是否存在未注释的节点
missing_nodes <- setdiff(V(gse20029_network)$name, gse20029_node_attributes$protein)
if (length(missing_nodes) > 0) {
  cat("⚠️ 缺失注释的节点数:", length(missing_nodes), "\n")
} else {
  cat("✅ 所有节点均已注释。\n")
}

# 执行 PCA 分析（基于平移后的拓扑指标）
pca_matrix <- gse20029_results %>%
  dplyr::select(Degree_Shifted, Betweenness_Shifted, Eigenvector_Shifted, MCC_Shifted, EPC_Shifted) %>%
  as.matrix()
pca_gse20029 <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)

# 拓扑权重计算函数：利用 PCA loadings 和解释方差计算各指标权重
calculate_weights <- function(pca_result) {
  loading <- as.data.frame(pca_result$rotation)
  explained_var <- summary(pca_result)$importance[2, ]
  weights <- rowSums((loading^2) * explained_var)
  return(weights / sum(weights))
}

gse20029_weights <- calculate_weights(pca_gse20029)

# 组织权重数据框
weight_df <- data.frame(
  Method = c("Degree_Zscore", "Betweenness_Zscore", "Eigenvector_Zscore", "MCC_Zscore", "EPC_Zscore"),
  GSE20029_Weight = gse20029_weights
)
print(weight_df)

# 保存权重表
write.table(weight_df,
            file = file.path(work_path, "results/PPI/GSE20029_Topology_Weights.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 可视化权重
weights_long <- melt(weight_df, id.vars = "Method")
ggplot(weights_long, aes(x = Method, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "GSE20029 Topology Scoring Weights", y = "Weight", x = "Method", fill = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = "#A064C7")

# 保存图像（确保目标文件夹存在）
fig_output_path <- file.path(work_path, "figures/Drug_Mapping")
dir.create(fig_output_path, showWarnings = FALSE, recursive = TRUE)
ggsave(filename = file.path(fig_output_path, "Weight_GSE20029_Barplot.png"),
       width = 10, height = 6, dpi = 100)

# 保存分析对象
save(gse20029_network, gse20029_results, weight_df, pca_gse20029,
     file = file.path(work_path, "results/PPI/TopoNetwork_GSE20029.RData"))

cat("✅ GSE20029 网络拓扑分析完成！结果已保存。\n")
