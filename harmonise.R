###### 差异蛋白质分析######
# 加载 limma 包
library(limma)

#设置路径
setwd("D:/postgraduate.research/hzx’s subject/excel proteome")
getwd()
#导入
proteome_log_file_copy_REH <- read.csv("proteome_log_file_copy_REH.csv",header = TRUE, check.names = FALSE)
#去除含缺失值#NAME?的行,将填充值设为numeric
proteome_log_file_deleteNA_REH <- proteome_log_file_copy_REH[!apply(proteome_log_file_copy_REH == "#NAME?", 1, any), ]
write.csv(proteome_log_file_deleteNA_REH,file = "proteome_log_file_deleteNA_REH.csv")
proteome_log_file_deleteNA_REH[,3:6] <- as.data.frame(lapply(proteome_log_file_deleteNA_REH[,3:6], as.numeric))
is.numeric(proteome_log_file_deleteNA_REH[,4])

#创建分组数据框
GroupFrame <- data.frame(
  group = c("ctrl", "treat", "ctrl", "treat"),
  row.names = c("REH", "REH_SCFV", "REH_293T", "REH_SCFV_copy")
)
group_list <- factor(GroupFrame[,1],ordered = F)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(proteome_log_file_deleteNA_REH[,3:6])

# 线性模型拟合
fit <- lmFit(proteome_log_file_deleteNA_REH[,3:6],design)

# 构建比较矩阵（分析的时候最好输出一下构建的比较矩阵，反了的话有点麻烦，这里系数大于0，则treat表达量高于ctrl）
cont.matrix <- makeContrasts(contrasts = paste0(rev(unique(group_list)),collapse = "-"),levels = design)
cont.matrix
fit2 <- contrasts.fit(fit,cont.matrix) # 构建数据线性模型,计算估计的相关系数和标准差
fit2 <- eBayes(fit2) # 贝叶斯检验
tmpOut <- topTable(fit2,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC") # 结果
limma.na <- na.omit(tmpOut) # 去除带有NA值的结果

#筛选差异蛋白
dif1 <- limma.na[limma.na$P.Value <= 0.05 & abs(limma.na$logFC) > log2(2),]
dif2 <- limma.na[limma.na$adj.P.Val <= 0.05 & abs(limma.na$logFC) > log2(2),]

#输出文件
write.csv(dif1,file = "dif1_REH_protein.csv")
write.csv(dif2,file = "dif2_REH_protein.csv")

#去excel中增加id列名，根据id merge
proteome_log_file_deleteNA_REH<-read.csv("proteome_log_file_deleteNA_REH.csv",header = TRUE, check.names = FALSE)
dif1<-read.csv("dif1_REH_protein.csv",header = TRUE, check.names = FALSE)
dif2<-read.csv("dif2_REH_protein.csv",header = TRUE, check.names = FALSE)
proteome_log_file_deleteNA_REH_dif1<-merge(proteome_log_file_deleteNA_REH,dif1,by = "id" ,all = F)
proteome_log_file_deleteNA_REH_dif2<-merge(proteome_log_file_deleteNA_REH,dif2,by = "id" ,all = F)
write.csv(proteome_log_file_deleteNA_REH_dif1,file="proteome_log_file_deleteNA_REH_dif1.csv")
write.csv(proteome_log_file_deleteNA_REH_dif2,file="proteome_log_file_deleteNA_REH_dif2.csv")

# 使用 ifelse() 设置 change 列的值
proteome_log_file_deleteNA_REH_dif1 <- read.csv("proteome_log_file_deleteNA_REH_dif1.csv",header = TRUE, check.names = FALSE)
proteome_log_file_deleteNA_REH_dif1$change <- ifelse(proteome_log_file_deleteNA_REH_dif1$t > 0, "up", "down")
write.csv(proteome_log_file_deleteNA_REH_dif1,file="proteome_log_file_deleteNA_REH_dif1.csv")
proteome_log_file_deleteNA_REH_dif1_up <- subset(proteome_log_file_deleteNA_REH_dif1, change == "up")
write.csv(proteome_log_file_deleteNA_REH_dif1_up,file="proteome_log_file_deleteNA_REH_dif1_up.csv")
proteome_log_file_deleteNA_REH_dif1_down <- subset(proteome_log_file_deleteNA_REH_dif1, change == "down")
write.csv(proteome_log_file_deleteNA_REH_dif1_down,file="proteome_log_file_deleteNA_REH_dif1_down.csv")