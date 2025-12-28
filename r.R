# ==============================================================================
# GO Biological Process Enrichment Visualization Script
# ==============================================================================
# 
# Description: This script visualizes DAVID GO enrichment analysis results for
#              hippocampal tissue transcriptomics study (LFD vs HFD mice)
# 
# Input: CSV file from DAVID GO enrichment analysis
# Output: Publication-quality figures (PNG, PDF, PPT)
# 
# Dependencies: readr, ggplot2, dplyr, tidyr, stringr, ggprism, ggthemes, forcats
# Optional: officer, rvg (for PPT generation)
# 
# Author: [Your Name]
# Date: [Date]
# Version: 1.0
# 
# Citation: If used in your research, please cite relevant tools:
#           DAVID (Huang et al., 2009), ggplot2 (Wickham, 2016)
# 
# Repository: https://github.com/[username]/[repository]
# 
# ==============================================================================

# 检查并安装必要的包
required_packages <- c("readr", "ggplot2", "dplyr", "tidyr", "stringr", 
                       "ggprism", "ggthemes", "forcats")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
}

# 安装缺失的包
invisible(sapply(required_packages, install_if_missing))

# 加载包
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggprism)
library(ggthemes)
library(forcats)

# 可选：加载PPT生成包（如果可用）
ppt_available <- requireNamespace("officer", quietly = TRUE) && 
  requireNamespace("rvg", quietly = TRUE)

if (ppt_available) {
  library(officer)
  library(rvg)
  cat("PPT generation packages loaded successfully.\n")
} else {
  cat("PPT generation packages not available. Skipping PPT creation.\n")
}

# ==============================================================================
# 参数设置
# ==============================================================================

# 输入文件路径（修改为您的文件路径）
csv_file <- "C:/Users/summe/OneDrive/Desktop/工作簿 2.csv"

# 输出目录
output_dir <- "C:/Users/summe/OneDrive/Desktop/"

# 预设的GO术语列表（23个术语）
specified_terms <- c(
  # Stress and Chemical Response (7个)
  "negative regulation of response to stimulus",
  "response to oxygen-containing compound",
  "cellular response to chemical stimulus",
  "multicellular organismal response to stress",
  "response to lipid",
  "cellular response to xenobiotic stimulus",
  "response to fatty acid",
  
  # Cell Proliferation and Death Regulation (4个)
  "regulation of neuron apoptotic process",
  "regulation of apoptotic process",
  "negative regulation of cell growth",
  "regulation of extrinsic apoptotic signaling pathway in absence of ligand",
  
  # Adhesion and Cytoskeleton Organization (7个)
  "regulation of chemotaxis",
  "regulation of cell adhesion",
  "regulation of cell-substrate adhesion",
  "regulation of cell-matrix adhesion",
  "regulation of cell-cell adhesion",
  "regulation of actin cytoskeleton organization",
  "positive regulation of extracellular matrix organization",
  
  # Cell Signal Transduction (5个)
  "intracellular signal transduction",
  "regulation of MAPK cascade",
  "regulation of ERK1 and ERK2 cascade",
  "regulation of protein kinase activity",
  "regulation of G protein-coupled receptor signaling pathway"
)

# 颜色方案
pal <- c(
  "Stress and\nChemical Response" = "#E41A1C",     # 红色
  "Cell Proliferation\nand Death Regulation" = "#377EB8", # 蓝色
  "Adhesion and\nCytoskeleton Organization" = "#4DAF4A",  # 绿色
  "Cell Signal\nTransduction" = "#984EA3"           # 紫色
)

# ==============================================================================
# 数据加载和预处理
# ==============================================================================

cat("===========================================\n")
cat("GO Enrichment Visualization Pipeline\n")
cat("===========================================\n\n")

# 1. 检查文件是否存在
if (!file.exists(csv_file)) {
  cat("错误: 文件不存在！请检查路径。\n")
  cat("文件路径:", csv_file, "\n")
  stop("无法找到文件")
}

cat("正在读取文件:", csv_file, "\n")

# 2. 读取CSV文件
go_bp_data <- read.table(
  csv_file,
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE,
  encoding = "UTF-8",
  fill = TRUE,
  comment.char = "",
  strip.white = TRUE
)

# 显示原始数据信息
cat("原始数据维度:", dim(go_bp_data), "行 x 列\n")
cat("原始列名:", paste(colnames(go_bp_data), collapse = ", "), "\n")

# 3. 简单处理数据 - 不删除任何行
go_bp_data <- go_bp_data[, 1:4]
colnames(go_bp_data) <- c("Term", "Count", "PValue", "Genes")

cat("\n处理后的数据维度:", dim(go_bp_data), "行 x 列\n")

# 4. 修复问号问题
go_bp_data$Term <- gsub("[?]", " ", go_bp_data$Term)
go_bp_data$Genes <- gsub("[?]", ", ", go_bp_data$Genes)
go_bp_data$Term <- trimws(go_bp_data$Term)
go_bp_data$Genes <- trimws(go_bp_data$Genes)

# 5. 数据转换
cat("\n开始数据转换...\n")

# 转换数值列
go_bp_data <- go_bp_data %>%
  mutate(
    Count = suppressWarnings(as.numeric(Count)),
    Count = ifelse(is.na(Count), 1, Count),
    
    PValue = suppressWarnings(as.numeric(PValue)),
    PValue = ifelse(is.na(PValue) | PValue <= 0, 0.5, PValue),
    
    neg_log10_pval = -log10(PValue)
  )

cat("转换后数据行数:", nrow(go_bp_data), "行\n")

# ==============================================================================
# GO术语筛选和分类
# ==============================================================================

cat("\n开始精确筛选和分类GO条目...\n")

# 首先清洗Term列，去掉GO编号部分
go_bp_data <- go_bp_data %>%
  mutate(
    Term_clean = gsub("^GO:\\d+~", "", Term),
    Term_clean = trimws(Term_clean)
  )

# 只保留指定的term
go_bp_data_filtered <- go_bp_data %>%
  filter(Term_clean %in% specified_terms)

cat("筛选后数据行数:", nrow(go_bp_data_filtered), "行\n")

# 验证是否所有23个术语都存在
if (nrow(go_bp_data_filtered) != 23) {
  cat("警告: 筛选后只有", nrow(go_bp_data_filtered), "个term，不是23个\n")
  cat("请检查Term名称是否完全匹配\n")
  
  # 找出缺失的term
  missing_terms <- setdiff(specified_terms, go_bp_data_filtered$Term_clean)
  if (length(missing_terms) > 0) {
    cat("缺失的term:\n")
    for (term in missing_terms) {
      cat("  -", term, "\n")
    }
  }
  
  cat("\n继续处理现有数据...\n")
}

# 根据term内容进行分类
go_bp_data <- go_bp_data_filtered %>%
  mutate(
    Subcategory = case_when(
      # Stress and Chemical Response (7个)
      Term_clean %in% c(
        "negative regulation of response to stimulus",
        "response to oxygen-containing compound",
        "cellular response to chemical stimulus",
        "multicellular organismal response to stress",
        "response to lipid",
        "cellular response to xenobiotic stimulus",
        "response to fatty acid"
      ) ~ "Stress and\nChemical Response",
      
      # Cell Proliferation and Death Regulation (4个)
      Term_clean %in% c(
        "regulation of neuron apoptotic process",
        "regulation of apoptotic process",
        "negative regulation of cell growth",
        "regulation of extrinsic apoptotic signaling pathway in absence of ligand"
      ) ~ "Cell Proliferation\nand Death Regulation",
      
      # Adhesion and Cytoskeleton Organization (7个)
      Term_clean %in% c(
        "regulation of chemotaxis",
        "regulation of cell adhesion",
        "regulation of cell-substrate adhesion",
        "regulation of cell-matrix adhesion",
        "regulation of cell-cell adhesion",
        "regulation of actin cytoskeleton organization",
        "positive regulation of extracellular matrix organization"
      ) ~ "Adhesion and\nCytoskeleton Organization",
      
      # Cell Signal Transduction (5个)
      Term_clean %in% c(
        "intracellular signal transduction",
        "regulation of MAPK cascade",
        "regulation of ERK1 and ERK2 cascade",
        "regulation of protein kinase activity",
        "regulation of G protein-coupled receptor signaling pathway"
      ) ~ "Cell Signal\nTransduction",
      
      TRUE ~ "Other"
    )
  )

# 显示分类结果
cat("\n分类结果:\n")
category_counts <- table(go_bp_data$Subcategory)
print(category_counts)

# ==============================================================================
# 数据排序和准备
# ==============================================================================

cat("\n开始排序...\n")

# 首先确保Subcategory是因子，按照需要的顺序
go_bp_data$Subcategory <- factor(go_bp_data$Subcategory, 
                                 levels = c("Stress and\nChemical Response", 
                                            "Cell Proliferation\nand Death Regulation", 
                                            "Adhesion and\nCytoskeleton Organization", 
                                            "Cell Signal\nTransduction"))

# 在每个亚组内按P值排序
go_bp_data <- go_bp_data %>%
  group_by(Subcategory) %>%
  arrange(PValue, .by_group = TRUE) %>%
  mutate(
    rank_in_group = row_number()
  ) %>%
  ungroup() %>%
  # 按亚组顺序和组内排序排列
  arrange(Subcategory, rank_in_group) %>%
  mutate(
    # 创建显示顺序（从上到下）
    display_order = rev(row_number()),  # 反转顺序，确保从上到下
    # 为每个条目创建唯一的因子
    term_factor = factor(Term_clean, levels = unique(Term_clean[order(display_order)]))
  )

# 创建亚组标签数据
rect.data <- go_bp_data %>%
  group_by(Subcategory) %>%
  summarise(
    min_idx = min(as.numeric(term_factor)),
    max_idx = max(as.numeric(term_factor)),
    count = n(),
    ycenter = (min_idx + max_idx) / 2
  ) %>%
  ungroup() %>%
  mutate(
    ymin = min_idx - 0.45,
    ymax = max_idx + 0.45
  )

# ==============================================================================
# 图形创建
# ==============================================================================

cat("\n开始创建图形...\n")

# 计算合适的x轴范围
max_pval <- max(go_bp_data$neg_log10_pval, na.rm = TRUE)
x_limit <- max_pval * 1.3  # 增加右侧空间

# 定义各个元素的位置
subcategory_x <- -max_pval * 0.6  # 亚组标签位置
count_circle_x <- -max_pval * 0.4  # 基因数量圆圈位置
term_name_x <- -max_pval * 0.3     # Term名称位置
bar_start_x <- 0                  # 条形图起始位置
gene_x <- bar_start_x + 0.1       # 基因列表位置

# 计算左侧限制
left_limit <- -max_pval * 0.8

# 为圆圈大小创建示例值 - 用于图例
count_range <- range(go_bp_data$Count)
# 创建几个代表性的值
legend_counts <- unique(round(seq(count_range[1], count_range[2], length.out = 4)))

# 创建图形
p <- ggplot(go_bp_data) +
  
  # 1. 亚组背景矩形
  geom_rect(
    aes(xmin = left_limit, xmax = subcategory_x + 0.2, 
        ymin = ymin, ymax = ymax, fill = Subcategory),
    data = rect.data,
    alpha = 0.15,
    color = NA
  ) +
  
  # 2. 亚组标签
  geom_text(
    aes(x = subcategory_x, y = ycenter, 
        label = paste0(Subcategory, "\n(n=", count, ")")),
    data = rect.data,
    size = 3.2,
    fontface = "bold",
    color = "black",
    lineheight = 0.8,
    hjust = 0.5
  ) +
  
  # 3. 基因数量圆圈 - 调大1.5倍
  geom_point(
    aes(x = count_circle_x, y = term_factor, size = Count),
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 1.0  # 稍微增加边框粗细
  ) +
  
  # 4. 基因数量数字标签
  geom_text(
    aes(x = count_circle_x, y = term_factor, label = Count),
    size = 2.8,
    fontface = "bold"
  ) +
  
  # 5. Term名称
  geom_text(
    aes(x = term_name_x, y = term_factor, label = Term_clean),
    hjust = 0,
    size = 3.0,
    fontface = "bold",
    color = "black"
  ) +
  
  # 6. 主条形图
  geom_col(
    aes(x = neg_log10_pval, y = term_factor, fill = Subcategory),
    width = 0.7,
    alpha = 0.85
  ) +
  
  # 7. 基因列表 - 在条形图内部，白色字体
  geom_text(
    aes(
      x = gene_x,
      y = term_factor,
      label = sapply(strsplit(Genes, ", "), function(x) {
        genes <- x[!is.na(x) & x != ""]
        if (length(genes) > 5) {
          paste0(paste(genes[1:5], collapse = ", "), "...")
        } else {
          paste(genes, collapse = ", ")
        }
      })
    ),
    hjust = 0,
    size = 2.4,
    fontface = "italic",
    color = "white",
    alpha = 0.95
  ) +
  
  # 8. 美化设置
  scale_fill_manual(values = pal, name = "Subcategory") +
  scale_size_continuous(
    name = "Gene Count",
    range = c(2.2 * 1.5, 4.5 * 1.5),  # 调大1.5倍
    breaks = legend_counts,
    guide = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "right",
      order = 2
    )
  ) +
  scale_x_continuous(
    name = expression(-log[10]("P Value")),
    expand = expansion(mult = c(0.05, 0.05)),
    limits = c(left_limit, x_limit),
    breaks = seq(0, ceiling(max_pval), by = 1)
  ) +
  scale_y_discrete(
    expand = expansion(add = 0.15)
  ) +
  
  # 9. 标签
  labs(
    title = "GO Biological Process Enrichment Analysis",
    subtitle = paste("Hippocampal tissue: LFD vs HFD | Showing", 
                     nrow(go_bp_data), "significant terms"),
    caption = "Data source: DAVID GO enrichment analysis | Visualization: R ggplot2"
  ) +
  
  # 10. 主题设置
  theme_prism() +
  theme(
    # 隐藏Y轴
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    
    # X轴设置
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
    
    # 标题
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      color = "gray40",
      margin = margin(b = 12)
    ),
    plot.caption = element_text(
      size = 9,
      hjust = 0.5,
      color = "gray50",
      margin = margin(t = 10)
    ),
    
    # 网格线
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor.x = element_blank(),
    
    # 图例
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = unit(0.3, "cm"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.5, "cm"),
    
    # 背景
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    
    # 边距
    plot.margin = margin(15, 40, 15, 15)
  ) +
  
  # 11. 图例顺序
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      order = 1
    )
  )

# 显示图形
cat("\n显示图形...\n")
print(p)

# ==============================================================================
# 保存图形
# ==============================================================================

cat("\n保存图形文件...\n")

# 创建输出目录（如果不存在）
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# 自动计算合适的高度
plot_height <- max(10, nrow(go_bp_data) * 0.5)

# 保存为PNG
png_file <- file.path(output_dir, "GO_BP_Enrichment.png")
ggsave(
  png_file,
  plot = p,
  width = 19,
  height = plot_height,
  dpi = 300,
  bg = "white"
)
cat("PNG图形已保存到:", png_file, "\n")
cat("图形尺寸:", 19, "x", round(plot_height, 1), "英寸\n")

# 保存为PDF
pdf_file <- file.path(output_dir, "GO_BP_Enrichment.pdf")
ggsave(
  pdf_file,
  plot = p,
  width = 19,
  height = plot_height,
  bg = "white"
)
cat("PDF图形已保存到:", pdf_file, "\n")

# ==============================================================================
# 可选：生成PPT文件
# ==============================================================================

if (ppt_available) {
  cat("\n正在生成PPT文件...\n")
  tryCatch({
    # 创建新的PowerPoint文档
    ppt <- read_pptx()
    
    # 添加幻灯片
    ppt <- ppt %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = "GO Biological Process Enrichment Analysis", 
              location = ph_location_type(type = "title")) %>%
      ph_with(value = dml(ggobj = p, 
                          bg = "white", 
                          width = 15,
                          height = plot_height * 0.9),
              location = ph_location_type(type = "body"))
    
    # 保存PPT文件
    ppt_file <- file.path(output_dir, "GO_BP_Enrichment.pptx")
    print(ppt, target = ppt_file)
    cat("PPT文件已保存到:", ppt_file, "\n")
  }, error = function(e) {
    cat("生成PPT文件时出错:", e$message, "\n")
  })
} else {
  cat("\n跳过PPT文件生成（缺少必需包）\n")
}

# ==============================================================================
# 生成分析报告
# ==============================================================================

cat("\n===========================================================\n")
cat("分析报告\n")
cat("===========================================================\n")

# 显示分类统计
cat("\n最终分类统计:\n")
final_stats <- go_bp_data %>%
  group_by(Subcategory) %>%
  summarise(
    Terms = n(),
    Total_Genes = sum(Count),
    Mean_Genes = round(mean(Count), 1),
    Min_PValue = formatC(min(PValue), format = "e", digits = 2),
    Max_PValue = formatC(max(PValue), format = "e", digits = 2),
    Mean_neg_log10_P = round(mean(neg_log10_pval), 2)
  ) %>%
  arrange(desc(Terms))

print(as.data.frame(final_stats))

# 总体统计
cat("\n总体统计:\n")
cat("总GO术语数:", nrow(go_bp_data), "\n")
cat("总基因数:", sum(go_bp_data$Count), "\n")
cat("平均每个术语基因数:", round(mean(go_bp_data$Count), 1), "\n")
cat("最小P值:", formatC(min(go_bp_data$PValue), format = "e", digits = 2), "\n")
cat("最大-log10(P值):", round(max(go_bp_data$neg_log10_pval), 2), "\n")

cat("\n===========================================================\n")
cat("分析完成！\n")
cat("===========================================================\n")

# 保存会话信息
session_file <- file.path(output_dir, "session_info.txt")
sink(session_file)
cat("GO Enrichment Analysis - Session Information\n")
cat("============================================\n")
cat("Analysis date:", date(), "\n\n")
cat("R version:\n")
print(R.version)
cat("\n\nLoaded packages:\n")
print(sessionInfo()$otherPkgs)
sink()

cat("会话信息已保存到:", session_file, "\n")

# ==============================================================================
# 文件列表
# ==============================================================================

cat("\n生成的文件列表:\n")
cat("1. PNG图像:", png_file, "\n")
cat("2. PDF文档:", pdf_file, "\n")
if (ppt_available && file.exists(file.path(output_dir, "GO_BP_Enrichment.pptx"))) {
  cat("3. PPT演示文稿:", file.path(output_dir, "GO_BP_Enrichment.pptx"), "\n")
}
cat("4. 会话信息:", session_file, "\n")

cat("\n===========================================================\n")
cat("脚本执行完成！\n")
cat("===========================================================\n")