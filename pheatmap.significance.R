library(pheatmap)
library(grid)
library(readxl)
install.packages("readr")
library(readr)
setwd("D:/postgraduate.research/my subject/Mendelian randomization/figure/pheatmap")

#读取数据
dat <- read.csv("pheatmap.csv",row.names = 1, header = T, check.names = F )
# 转置数据框
dat_t <- as.data.frame(t(dat))

# 保留原来的列名作为新的行名
colnames(dat_t) <- row.names(dat)

significance <- read.csv("significance.csv",row.names = 1, header = T, check.names = F )
significance_t <- as.data.frame(t(significance))
colnames(significance_t) <- row.names(significance)

# 对数转换
# +0.1是为了防止对0取对数；是加0.1还是加个更小的值取决于数据的分布。
# 加的值一般认为是检测的低阈值，低于这个值的数字之间的差异可以忽略。
#dat1_log <- log2(dat1+0.1)

#添加样本分组
#annotation_col1 = data.frame(
#  Group = c('control', 'test'),
#  row.names = colnames(dat1_log)
#)

#绘制基因表达热图
#以下仅列举简单参数，关于颜色、字体等的详细参数等参阅帮助 
#?pheatmap
# 输出R默认大小，分辨率增加到1000（默认69）

morandi_blue_red <- colorRampPalette(c('#3A6DAB', 'white', 'red'))(1000)
dev.new()
pheatmap(
  dat_t,
  cluster_rows = TRUE, cluster_cols = T,  #行列是否聚类
  fontsize_row = 7, fontsize_col = 7,  #字体大小设置
  color = morandi_blue_red,
  #color = colorRampPalette(c('blue', 'white', 'red'))(1000),  #热图颜色设置，本示例由绿到红渐变表示表达值增加
  #annotation_col = annotation_col1,  #定义样本分组
  cellwidth = 10, cellheight = 15, #border_color = "black", #格子宽度高度等设置
  show_rownames = TRUE,  # 显示行名 # 设置行名为 gene_symbol 列的值
  annotation_names_col = F,
  display_numbers = significance_t,
  #fontface="italic",
  fontfamily= "Times New Roman"
  #annotation_col = NULL,
  #filename ="DEG.DESeq2.up.all.REH.SEM.png"
  
)  #字体大小
#关闭当前图形设备
#dev.off()






