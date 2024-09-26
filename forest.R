#install.packages("forestplot")
library(forestplot)
library(readxl)
setwd("D:/postgraduate.research/my subject/Mendelian randomization/figure/forest")
options(scipen = 1)
options(digits = 2)
forest <- read_excel("CIF to HM forest map.xlsx")
row_name<-cbind(c("exposure",forest$V1),c("outcome",forest$V2),c("OR(95%CI)",forest$V3),c("P-value",forest$V7)) 
forest<-rbind(rep(NA,7),forest) 

#pdf(file="森林图.pdf",
#width = 10,             #图片的宽度
#   height =7,            #图片的高度
#)
#?forestplot
is.summary <- c(TRUE, rep(FALSE, nrow(forest) - 1))
dev.new()
forestplot(labeltext=row_name,forest[,c("V4","V5","V6")],
           zero = 1,xticks = c(0.1,1.0,2.5,10),boxsize = 0.1,lineheight = unit(6,"mm"),
           colgap = unit(.2,"npc"),lwd.zero = 0.5,lwd.ci = 2,lty.ci = 2,
           col = fpColors(box = "blue",summary = "black",zero = "red"),
           lwd.xaxis = 1,graph.pos = 3,graphwidth = unit(.4,"npc"),is.summary = is.summary,
           #xlog = TRUE, # 使x轴呈对数尺度
           clip = c(0.1, 5), # 将x轴的显示范围限制在0.1到10之间
          # txt_gp = fpTxtGp(label = gpar(fontsize = 12), ticks = gpar(fontsize = 10), xlab = gpar(fontsize = 14))
)# 图形文本部分 T所对应行加粗 graphwidth调整森林图所在表格的宽度,colgap = unit(.2,"npc")调整列间距
dev.off()

