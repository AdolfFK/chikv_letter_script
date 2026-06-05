#!/usr/bin/env Rscript

# ===============================
# 加载必要的包
# ===============================
suppressPackageStartupMessages({
    library(survival)
    library(survminer)
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(dplyr)
    library(grid)
    library(gridExtra)
    library(showtext)
})

# ===============================
# 字体设置
# ===============================
font_add("Arial", regular = "Arial.ttf")
showtext_auto()

# ===============================
# 路径设置
# ===============================
setwd("/Users/adolf1/Documents/Work/Data_wangxx/Chikungunya-yuecan/pipiline_v2/plot-20260206/4")
input_file <- "convert_data_v2/cox_ready_data.csv"
outdir <- "cox_analysis_v2_forest"
fig_dir <- file.path(outdir, "figures")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ===============================
# 读取并准备数据
# ===============================
cat("读取数据...\n")
df_raw <- read_csv(
  input_file,
  na = c("", "NA", "N/A", "None", "none"),
  show_col_types = FALSE
)

df <- df_raw %>%
  mutate(
    Age = as.numeric(Age),
    CT_Value = as.numeric(CT_Value),
    Illness_Duration_Days = as.numeric(Illness_Duration_Days),
    event = as.numeric(event),
    Gender = factor(Gender, levels = c("0", "1"), labels = c("Female", "Male")),
    Has_Fever = factor(Has_Fever, levels = c("0", "1"), labels = c("No", "Yes")),
    Has_Rash = factor(Has_Rash, levels = c("0", "1"), labels = c("No", "Yes")),
    Has_Joint_Pain = factor(Has_Joint_Pain, levels = c("0", "1"), labels = c("No", "Yes"))
  )

df_cox <- df %>%
  select(Age, Gender, CT_Value, Has_Fever, Has_Rash, Has_Joint_Pain, 
         Illness_Duration_Days, event) %>%
  na.omit()

surv_obj <- with(df_cox, Surv(Illness_Duration_Days, event))

# ===============================
# 多因素Cox回归分析
# ===============================
cat("拟合多因素Cox回归模型...\n")
multi_fit <- coxph(
  surv_obj ~ Gender + Age + CT_Value + Has_Fever + Has_Rash + Has_Joint_Pain,
  data = df_cox
)

# ===============================
# 导出统计表（核心新增功能）
# ===============================
cat("导出Cox结果表...\n")
model_summary <- summary(multi_fit)
results_table <- data.frame(
  Variable = rownames(model_summary$coefficients),
  HR = model_summary$coefficients[, "exp(coef)"],
  CI_Low = model_summary$conf.int[, "lower .95"],
  CI_High = model_summary$conf.int[, "upper .95"],
  P_value = model_summary$coefficients[, "Pr(>|z|)"]
)
# 排序 & 保持和图一致
variable_order <- c("GenderMale", "Age", "CT_Value", 
                    "Has_FeverYes", "Has_RashYes", "Has_Joint_PainYes")
results_table <- results_table %>%
  dplyr::filter(Variable %in% variable_order) %>%
  dplyr::mutate(Variable = factor(Variable, levels = variable_order)) %>%
  dplyr::arrange(Variable)
# 美化（可选）
results_table <- results_table %>%
  mutate(
    HR = round(HR, 3),
    CI_Low = round(CI_Low, 3),
    CI_High = round(CI_High, 3),
    P_value = signif(P_value, 3)
  )
# 输出目录
table_dir <- file.path(outdir, "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
# 保存 CSV
write.csv(results_table,
          file = file.path(table_dir, "cox_results_table.csv"),
          row.names = FALSE)
# 保存 TSV（论文常用）
write.table(results_table,
            file = file.path(table_dir, "cox_results_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("统计表已保存到: ", table_dir, "\n")

# ===============================
# 修改后的基础绘图函数 - 图例进一步向右移动
# ===============================
nature_forest_plot <- function(cox_model, df_cox) {
  
  model_summary <- summary(cox_model)
  
  results <- data.frame(
    Variable = rownames(model_summary$coefficients),
    HR = model_summary$coefficients[, "exp(coef)"],
    CI_low = model_summary$conf.int[, "lower .95"],
    CI_high = model_summary$conf.int[, "upper .95"],
    p_value = model_summary$coefficients[, "Pr(>|z|)"]
  )
  
  variable_order <- c("GenderMale", "Age", "CT_Value", 
                      "Has_FeverYes", "Has_RashYes", "Has_Joint_PainYes")
  
  results <- results %>%
    dplyr::filter(Variable %in% variable_order) %>%
    dplyr::mutate(Variable = factor(Variable, levels = rev(variable_order)))
  
  display_names <- c(
    "GenderMale" = "Gender (Male)",
    "Age" = "Age",
    "CT_Value" = "CT Value",
    "Has_FeverYes" = "Has Fever (Yes)",
    "Has_RashYes" = "Has Rash (Yes)",
    "Has_Joint_PainYes" = "Has Joint Pain (Yes)"
  )
  
  results$Label <- display_names[as.character(results$Variable)]
  
  results$p_text <- ifelse(results$p_value < 0.001, "<0.001",
                           sprintf("%.3f", results$p_value))
  
  results$stars <- dplyr::case_when(
    results$p_value < 0.001 ~ "***",
    results$p_value < 0.01 ~ "**",
    results$p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  results$HR_text <- sprintf("%.2f (%.2f–%.2f)%s",
                             results$HR, results$CI_low, results$CI_high,
                             results$stars)
  
  results$color <- ifelse(results$HR > 1, "#D55E00", "#0072B2")
  results$p_color <- ifelse(results$p_value < 0.05, "red", "black")
  
  # ===================== 森林图 =====================
  p_forest <- ggplot2::ggplot(results, ggplot2::aes(x = HR, y = Variable)) +
    
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                        color = "grey30", linewidth = 1.8) +
    
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = CI_low, xmax = CI_high, color = color),
      height = 0.25, linewidth = 2.2
    ) +
    
    ggplot2::geom_point(ggplot2::aes(color = color), size = 4.5) +
    
    ggplot2::scale_color_identity() +
    
    ggplot2::theme_classic(base_family = "Arial") +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      text = ggplot2::element_text(family = "Arial", size = 17),
      legend.position = "none"
    ) +
    
    ggplot2::xlab("Hazard Ratio (HR)")
  
  # ===================== 左侧变量（✅已修正对齐） =====================
  p_left <- ggplot2::ggplot(results, ggplot2::aes(y = Variable, x = 0.05)) +
    ggplot2::geom_text(
      ggplot2::aes(label = Label),
      hjust = 0,
      nudge_x = -0.02,
      size = 6
    ) +
    ggplot2::xlim(0, 1.2) +
    ggplot2::theme_void(base_family = "Arial")
  
  # ===================== 右侧表格 =====================
  p_right <- ggplot2::ggplot(results, ggplot2::aes(y = Variable)) +
    
    ggplot2::geom_text(
      ggplot2::aes(x = 0, label = p_text, color = p_color),
      hjust = 0, size = 6, show.legend = FALSE
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(x = 1.8, label = HR_text),
      hjust = 0, size = 6
    ) +
    
    ggplot2::scale_color_identity() +
    ggplot2::theme_void(base_family = "Arial") +
    ggplot2::xlim(0, 4.5)
  
  # ===================== 右侧综合图例 =====================
  p_legend <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    
    ggplot2::annotate("text", x = 0, y = 1, label = "Color Legend:",
                      hjust = 0, size = 7, family = "Arial", fontface = "bold") +
    
    ggplot2::annotate("point", x = 0, y = 0.9, color = "#D55E00", size = 5) +
    ggplot2::annotate("text", x = 0.15, y = 0.9, label = "HR > 1 (Risk)",
                      hjust = 0, size = 6, family = "Arial") +
    
    ggplot2::annotate("point", x = 0, y = 0.8, color = "#0072B2", size = 5) +
    ggplot2::annotate("text", x = 0.15, y = 0.8, label = "HR < 1 (Protective)",
                      hjust = 0, size = 6, family = "Arial") +
    
    ggplot2::annotate("text", x = 0, y = 0.65, label = "Significance:",
                      hjust = 0, size = 7, family = "Arial", fontface = "bold") +
    
    ggplot2::annotate("text", x = 0, y = 0.55, label = "*   p < 0.05",
                      hjust = 0, size = 6, family = "Arial") +
    ggplot2::annotate("text", x = 0, y = 0.48, label = "**  p < 0.01",
                      hjust = 0, size = 6, family = "Arial") +
    ggplot2::annotate("text", x = 0, y = 0.41, label = "*** p < 0.001",
                      hjust = 0, size = 6, family = "Arial") +
    
    ggplot2::annotate("text", x = 0, y = 0.25, label = "Sample Info:",
                      hjust = 0, size = 7, family = "Arial", fontface = "bold") +
    
    ggplot2::annotate("text", x = 0, y = 0.18,
                      label = paste0("Total n = ", nrow(df_cox)),
                      hjust = 0, size = 6, family = "Arial") +
    
    ggplot2::annotate("text", x = 0, y = 0.11,
                      label = paste0("Events = ", sum(df_cox$event)),
                      hjust = 0, size = 6, family = "Arial") +
    
    ggplot2::annotate("text", x = 0, y = 0.04,
                      label = "Ref: Female, No symptoms",
                      hjust = 0, size = 6, family = "Arial") +
    
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
  
  # ===================== 标题（✅左对齐修正） =====================
  title_left <- grid::textGrob(
    "Variable",
    x = 0,
    just = "left",
    gp = grid::gpar(fontfamily="Arial", fontsize=22, fontface="bold")
  )
  
  title_mid <- grid::textGrob(" ", gp = grid::gpar(fontfamily="Arial"))
  
  title_right <- grid::textGrob(
    "P Value        HR (95% CI)",
    gp = grid::gpar(fontfamily="Arial", fontsize=22, fontface="bold")
  )
  
  # ===================== 合并 =====================
  gridExtra::grid.arrange(
    gridExtra::arrangeGrob(title_left, p_left, ncol=1, heights=c(0.12,1)),
    gridExtra::arrangeGrob(title_mid, p_forest, ncol=1, heights=c(0.12,1)),
    gridExtra::arrangeGrob(title_right, p_right, ncol=1, heights=c(0.12,1)),
    p_legend,
    ncol = 4,
    widths = c(3.2, 4, 4, 3)
  )
}






# 保存修改后的森林图（进一步增加宽度）
pdf(file.path(fig_dir, "forest_legend_far_right.pdf"), width = 14, height = 10)
nature_forest_plot(multi_fit, df_cox)
dev.off()


png(file.path(fig_dir, "forest_legend_far_right.png"), width = 14, height = 10, 
    units = "in", res = 300)
nature_forest_plot(multi_fit, df_cox)
dev.off()


