library(tidyverse)
library(glue)
library(tidygraph)
library(ggraph)
library(igraph)
source("gather_graph_node.R")
source("gather_graph_edge.R")
# 加载数据
setwd("D:/postgraduate.research/my subject/Mendelian randomization/figure/circle bubble/bad")
Promoter <- read.csv("bad.csv", header = TRUE)
Promoter$value <- rep(1)

# 构建颜色分类数据集
Promoter_index <- c("Cluster", "Proteins")
nodes <- gather_graph_node(Promoter, index = Promoter_index, value = "value", root = "cluster")
edges <- gather_graph_edge(Promoter, index = Promoter_index, root = "cluster")

# 构建图形对象
graph_cluster <- tbl_graph(nodes, edges)

# 创建圆形树状图
gc <- ggraph(graph_cluster, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(aes(color = node1.node.branch), alpha = 1/3) + 
  geom_node_point(aes(size = node.size, color = node.branch), alpha = 1/3) + 
  coord_fixed() +
  theme_void() + # 移除背景网格和轴线
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))

# 外环文字
gc <- gc + geom_node_text(
  aes(
    x = 1.02 * x,
    y = 1.02 * y,
    label = node.short_name,
    angle = -((-node_angle(x, y) + 90) %% 180) + 90,
    filter = leaf,
    color = node.branch
  ),
  size = 3, hjust = 'outward'
)

# 添加中心文本
gc <- gc + 
  annotate("text", x = 0, y = 0, label = "Promoter", size = 3.5, fontface = "bold", color = "#B22222")

# 调整主题设置，移除黑色框框
gc + theme_void() + theme(legend.position = "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())