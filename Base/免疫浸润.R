###

###
###

###

### install packages for analysis
depens<-c('tibble', 'survival', 'survminer', 
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', "biomaRt",
          'ggpubr',  "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}
if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR")

### load the packages for analysis -------------
library(IOBR)

count <- read.csv("C:/Users/ZHH/Desktop/G12D/lung.csv",row.names = 1)

expr_coad <- log2(count+0.1)
expr_coad <- expr_coad[apply(expr_coad,1,sd)>0.5,]

# MCPcounter
im_mcpcounter <- deconvo_tme(eset = expr_coad,
                             method = "mcpcounter"
)
# EPIC
im_epic <- deconvo_tme(eset = expr_coad,
                       method = "epic",
                       arrays = F
)
# xCell
im_xcell <- deconvo_tme(eset = expr_coad,
                        method = "xcell",
                        arrays = F
)
# CIBERSORT
im_cibersort <- deconvo_tme(eset = expr_coad,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)
# IPS
im_ips <- deconvo_tme(eset = expr_coad,
                      method = "ips",
                      plot = F
)
## 
## >>> Running Immunophenoscore

# quanTIseq
im_quantiseq <- deconvo_tme(eset = expr_coad,
                            method = "quantiseq",
                            scale_mrna = T
)
# ESTIMATE
im_estimate <- deconvo_tme(eset = expr_coad,
                           method = "estimate"
)
## 
## >>> Running ESTIMATE
# TIMER
im_timer <- deconvo_tme(eset = expr_coad
                        ,method = "timer"
                        ,group_list = rep("coad",dim(expr_coad)[2])
)
library(tidyr)
# 取前12个样本做演示
res<-cell_bar_plot(input = im_cibersort, title = "CIBERSORT Cell Fraction")
tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by="ID") %>% 
  inner_join(im_xcell, by="ID") %>% 
  inner_join(im_cibersort, by="ID") %>% 
  inner_join(im_ips, by= "ID") %>% 
  inner_join(im_quantiseq, by="ID") %>% 
  inner_join(im_estimate, by= "ID") %>% 
  inner_join(im_timer, by= "ID")