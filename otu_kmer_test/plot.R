library(vegan)
library(ggplot2)
tab <- read.table("/Users/wangpeng/Desktop/out.mat", comment.char = "",
                  header = T, sep = "\t", row.names = 1)

map <- read.table("/Users/wangpeng/Desktop/otu.map", comment.char = "")

tab2 <- t(tab)

pca <- rda(tab2)
d <- pca$CA$u
pdt <- data.frame(pc1=d[,1], pc2=d[,2],
                  group=map$V2[match(rownames(d), map$V1)])
rownames(pdt) <- rownames(pca$li)

p <- ggplot(pdt, aes(pc1, pc2, group = group, fill = group, color = group)) + 
  geom_point() 

ggsave(p, filename = "/Users/wangpeng/Desktop/pca_tnf.png")


