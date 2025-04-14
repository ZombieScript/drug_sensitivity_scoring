rm= ls()

library(readxl)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(tibble)

# ==== 参数 ====
n <- 12
obs_successes <- 8
p_null <- 0.5
p_alt <- obs_successes / n  # 观测准确率作为 H1

# ==== x 值 ====
x_vals <- 0:n
h0_probs <- dbinom(x_vals, n, p_null)   # H0 分布（蓝）
h1_probs <- dbinom(x_vals, n, p_alt)    # H1 分布（红）

# ==== 合并分布数据 ====
h1_label <- paste0("H1: p = ", round(p_alt, 3))  # label 用作分组
df_all <- tibble(
  hits = rep(x_vals, 2),
  probability = c(h0_probs, h1_probs),
  Distribution = rep(c("H0: p = 0.5", h1_label), each = length(x_vals))
)

# ==== P 值区域：H0 下 >= observed 的概率 ====
sig_region_x <- x_vals[x_vals >= obs_successes]
df_sig <- tibble(
  hits = sig_region_x,
  prob = dbinom(sig_region_x, n, p_null)
)

# ==== Power 区域：H1 分布落入显著区域的概率 ====
df_power <- tibble(
  hits = sig_region_x,
  prob = dbinom(sig_region_x, n, p_alt)
)

# ==== 计算数值 ====
pval <- signif(sum(df_sig$prob), 3)
power_val <- signif(sum(df_power$prob), 3)

# ==== 自定义颜色映射 ====
color_vals <- setNames(c("blue", "red"), c("H0: p = 0.5", h1_label))

# ==== 绘图 ====
p <- ggplot() +
  # 1. 灰色区域：H0 下 ≥ obs 的概率（P 值区域）
  geom_col(data = df_sig, aes(x = hits, y = prob), fill = "gray80", width = 0.8) +
  
  # 2. 红色透明区域：H1 下命中显著区的概率（Power 区域）
  geom_col(data = df_power, aes(x = hits, y = prob), fill = "red", alpha = 0.3, width = 0.8) +
  
  # 3. 两条分布曲线
  geom_line(data = df_all, aes(x = hits, y = probability, color = Distribution, linetype = Distribution), size = 1) +
  geom_point(data = df_all, aes(x = hits, y = probability, color = Distribution), size = 2) +
  
  # 4. 虚线标注：H0 期望 & 观测值
  geom_vline(xintercept = n * p_null, linetype = "dashed", color = "blue", size = 1) +
  geom_vline(xintercept = obs_successes, linetype = "dashed", color = "red", size = 1) +
  
  # 5. 注释 P 值 和 Power
  annotate("text", x = n - 0.3, y = max(h0_probs) * 0.9, 
           label = paste0("p = ", pval, "\nPower = ", power_val),
           hjust = 1, size = 5) +
  
  # 6. 坐标轴 & 样式
  scale_x_continuous(breaks = 0:n) +
  scale_color_manual(values = color_vals) +
  labs(
    title = "Binomial Test with Significance and Power Regions",
    subtitle = paste0("Sample size = ", n, ", Observed Hits = ", obs_successes),
    x = "Number of Correct Predictions",
    y = "Probability"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

# ==== 保存图像 ====
ggsave("E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/figures/Binomial_with_Power_and_Significance.png",
       p, width = 9, height = 5.5, dpi = 300)

# ==== 1. 读取数据 ====
file_path <- "E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/Summary of Datasets accuracy.csv"
df <- read_excel(file_path)

# ==== 2. 数据处理 ====
df <- df %>%
  mutate(
    Relative_Ranking_Pct = as.numeric(`Relative_Ranking(%)`),
    Predictive_Accuracy = Predictive_Accuracy == "Yes",
    Hit = Predictive_Accuracy & Relative_Ranking_Pct <= 10,
    Label = paste(Drugs, "|", Actual_Effect, "|", Cell_Lines)
  ) %>%
  arrange(desc(!is.na(Relative_Ranking_Pct)), Relative_Ranking_Pct) %>%  # 保证 NA 在最后
  mutate(Label = factor(Label, levels = Label))  # 保持顺序

# ==== 3. 二项分布检验 ====
total <- nrow(df)  # 包括 NA
successes <- sum(df$Hit, na.rm = TRUE)
binom_res <- binom.test(successes, total, p = 0.1, alternative = "greater")
pval <- binom_res$p.value
p_text <- ifelse(pval < 0.05, paste0("p = ", signif(pval, 3), " *"), paste0("p = ", signif(pval, 3)))

# # ==== 设定 P 值在右下角位置 ====
# ==== 设定 P 值显示在右下角 ====
# 注意：由于 coord_flip 翻转了，所以要把“第一个 label”当作最下方的那一条
x_pos <- df$Label[1]  # 图的最下方（左侧 y 轴）
y_pos <- max(df$Relative_Ranking_Pct, na.rm = TRUE) * 0.95  # 靠近右侧

# ==== 绘图 ====
p2 <- ggplot(df, aes(x = Label, y = Relative_Ranking_Pct, fill = Hit)) +
  geom_col(color = "black", na.rm = TRUE) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(
    title = "Drug Ranking and Prediction Accuracy",
    x = NULL,
    y = "Relative Ranking (%)"
  ) +
  annotate("text", x = x_pos, y = y_pos, label = p_text, hjust = 1, size = 5) +  # ✅ 修正位置
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.title = element_blank()
  )


# ==== 5. 绘图 ====
p2 <- ggplot(df, aes(x = Label, y = Relative_Ranking_Pct, fill = Hit)) +
  geom_col(color = "black", na.rm = TRUE) +  # 忽略 NA 高度
  geom_hline(yintercept = 10, linetype = "dashed", color = "blue") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(
    title = "Drug Ranking and Prediction Accuracy",
    x = NULL,
    y = "Relative Ranking (%)"
  ) +
  annotate("text", x = x_pos, y = y_pos, label = p_text, hjust = 1, size = 5) +  # 放右下角
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.title = element_blank()
  )

# ==== 6. 保存图像 ====
ggsave("E:/Bioinformatics/CDK4&6 resistance 2cell lines/Method_Validation/figures/Figure_2A_Final.png",
       p2, width = 10, height = 6, dpi = 300)