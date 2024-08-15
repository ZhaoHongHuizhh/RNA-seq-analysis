library(TCGAbiolinks)
library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)

# 下载并准备突变数据
mut_query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(mut_query)
mut_data <- GDCprepare(mut_query)

# 下载并准备基因表达数据
exp_query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(exp_query)
exp_data <- GDCprepare(exp_query)

# 提取ARID1A突变信息
arid1a_mut_samples <- mut_data %>%
  filter(Hugo_Symbol == "ARID1A") %>%
  select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  pull(Tumor_Sample_Barcode)

# 将样本分组
mut_samples <- unique(mut_data$Tumor_Sample_Barcode)
arid1a_status <- ifelse(mut_samples %in% arid1a_mut_samples, "Mutant", "Wildtype")

# 创建样本分组数据框
sample_groups <- data.frame(
  Sample = mut_samples,
  ARID1A_Status = arid1a_status
)

# 查看 exp_data 对象中可用的 assay 名称
assay_names <- assayNames(exp_data)
print(assay_names)

# 假设 assay 名称是 "htseq_counts"
# 提取表达数据矩阵
exp_matrix <- assay(exp_data, "htseq_counts")

# 检查表达矩阵的前几行和前几列
head(exp_matrix)

# 确保样本ID匹配
sample_groups <- sample_groups %>%
  filter(Sample %in% colnames(exp_matrix))

# 重新排序表达矩阵以匹配分组信息
exp_matrix <- exp_matrix[, sample_groups$Sample]

# 创建DESeq2数据对象
colData <- data.frame(
  row.names = sample_groups$Sample,
  ARID1A_Status = factor(sample_groups$ARID1A_Status)
)
dds <- DESeqDataSetFromMatrix(countData = exp_matrix, colData = colData, design = ~ ARID1A_Status)

# 运行差异表达分析
dds <- DESeq(dds)
res <- results(dds)

# 提取显著差异表达基因
sig_genes <- res %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  arrange(padj)

# 将显著差异表达基因转换为ENTREZID
sig_genes_entrez <- bitr(rownames(sig_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# 进行GO富集分析
go_enrichment <- enrichGO(
  gene = sig_genes_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# 进行KEGG富集分析
kegg_enrichment <- enrichKEGG(
  gene = sig_genes_entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# 输出结果
head(go_enrichment)
head(kegg_enrichment)
