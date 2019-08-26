# R script to analyze bicmix output

source('../SSLB_functions.R', echo=TRUE)

dir <- "bicmix_result"

library(pals)

N <- 3005
G <- 5000

meta_data <- read.table(file = "data/meta_data_mRNA.txt", stringsAsFactors = F)

cells <- read.table(file = "data/cell_info.txt", stringsAsFactors = F)
genes <- read.table(file = "data/gene_info.txt", stringsAsFactors = F)

cells <- as.matrix(cells)
genes <- as.matrix(genes)

K_init <- 100

X_bicmix <- try(as.matrix(read.table(paste(dir, "/EX_10000", sep = ""), header = F)), silent = TRUE)
if (inherits(X_bicmix, 'try-error')) {
  X_bicmix <- matrix(0, nrow = K_init, ncol = N)
}
X_bicmix <- t(X_bicmix)
B_bicmix <- try(as.matrix(read.table(paste(dir, "/LAM_10000", sep = ""), header = F)), silent = TRUE)
if (inherits(B_bicmix, 'try-error')) {
  B_bicmix <- matrix(0, nrow = G, ncol = K_init)
}

if (any(X_bicmix != 0)) {
  X_norm <- apply(X_bicmix, 2, function(x) sum(abs(x)))
  B_norm <- apply(B_bicmix, 2, function(x) sum(abs(x)))
  
  d <- sqrt(X_norm/B_norm)
  X_bicmix <- t(1/d * t(X_bicmix))
  B_bicmix <- t(d * t(B_bicmix))
  
  X_bicmix[abs(X_bicmix) < 10^(-10)] <- 0
  
}

K_bicmix <- ncol(X_bicmix)

# PLOT

X_bicmix[X_bicmix > 3] <- 3
X_bicmix[X_bicmix < -3] <- -3


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
  theme(panel.border = element_blank()) 

reorder <- c(union(c(which(apply(X_bicmix, 2, function(x) all(x!=0))),
             62, 64, 78,7,
             71, 72, 82, 91,
             55, 57, 
             65, 8, 32, 33, 87, 
             1, 10, 24, 35, 49, 11, 39,63, 66, 74, 9,
             89, 92, 44, 25, 13, 17, 67,
             69, 70
             ), (1:K_bicmix)[-73]), 73)
             

bicmix_plot <- plot_matrix(X_bicmix[ord, reorder], title = "BicMix Factor Matrix")


pdf(file = "figures/zeisel_bicmix.pdf", width = 10, height = 8, pointsize = 16)

grid.arrange(group_plot, subclass_plot, bicmix_plot, nrow = 1, widths = c(5, 1, 25))

dev.off()
