
source('../SSLB_functions.R', echo=TRUE)

library(pals)

load("ISA_result/out_ISA.RData")

X_ISA = out_ISA$rows
B_ISA = out_ISA$columns

K_ISA = ncol(X_ISA)

N <- 3005
G <- 5000

meta_data <- read.table(file = "data/meta_data_mRNA.txt", stringsAsFactors = F)

cells <- read.table(file = "data/cell_info.txt", stringsAsFactors = F)
genes <- read.table(file = "data/gene_info.txt", stringsAsFactors = F)

cells <- as.matrix(cells)
genes <- as.matrix(genes)

groups <- meta_data$group
subclasses <- as.numeric(factor(meta_data$level2class, levels = unique(meta_data$level2class)))
subclass_names <- as.character(meta_data$level2class)

group_names <- c("Interneurons", "S1 Pyramidal", "CA1 Pyramidal", 
                 "Oligodendrocytes", "Microglia", "Endothelial",
                 "Astrocytes", "Ependymal", "Mural")

# this subclass has no category
no_cat <- unique(subclasses[meta_data$level2class == "(none)"])

# order by group and subclass given group

ord <- order(subclasses[which(groups == 1)])
old <- which(groups == 1)

for (i in 2:max(groups)) {
  cur <- which(groups == i)
  add <- max(old) + order(subclasses[cur])
  old <- cur
  ord <- c(ord, add)
}

n_groups <- max(groups)
label_data <- data.frame(x = rep(-5, n_groups), 
                         y = c(2800, 2500, 1800, 1000, 520, 370, 200, 70, 10),
                         labels = group_names)

group_plot <- plot_matrix(groups[ord], title = "") + 
  scale_fill_gradientn(colors = brewer.paired(n_groups)) +
  theme(plot.margin=unit(c(5.5, 5.5, 5.5, 45),"pt")) +
  geom_text(data = label_data, aes(x = x, y = y, label = labels)) +
  coord_cartesian(clip = 'off') +
  theme(panel.border =element_blank())

subclass_cols <- numeric(max(subclasses))
subclass_cols[-c(no_cat)] <- brewer.paired(47)
subclass_cols[no_cat] <- "black"

subclass_plot <- plot_matrix(subclasses[ord], title = "") + 
  scale_fill_gradientn(colors = subclass_cols) +
  theme(panel.border =element_blank()) 


reorder = c()
for (i in 1:length(groups)) {
  reorder = c(reorder, which(apply(X_ISA[which(groups == i),], 2, function(x) sum(x!=0)) > 11))
}
reorder = union(reorder, 1:K_ISA)

pdf(file = "figures/zeisel_ISA.pdf", width = 10, height = 8, pointsize = 16)

grid.arrange(group_plot, subclass_plot, plot_matrix(X_ISA[ord,reorder], title = ""), nrow = 1, widths=c(5, 1, 25))

dev.off()




