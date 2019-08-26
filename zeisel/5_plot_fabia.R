# R script to analyze fabia biclusters

source('../SSLB_functions.R', echo=TRUE)

library(pals)

meta_data <- read.table(file = "data/meta_data_mRNA.txt", stringsAsFactors = F)

load("fabia_output.RData")

N <- 3005
G <- 5000

K_FABIA <- 100

X_FABIA <- t(out_FABIA@Z)
B_FABIA <- out_FABIA@L

FABIA_biclusters <- extractBic(out_FABIA)
B_FABIA_support <- matrix(0, nrow = G, ncol = K_FABIA)
X_FABIA_support <- matrix(0, nrow = N, ncol = K_FABIA)

for(k in 1:K_FABIA) {
  B_FABIA_support[FABIA_biclusters$numn[k,]$numng, k] <- 1
  X_FABIA_support[FABIA_biclusters$numn[k,]$numnp, k] <- 1
}

X_FABIA[X_FABIA_support == 0] <- 0
B_FABIA[B_FABIA_support == 0] <- 0

zeroes_B <- which(apply(B_FABIA, 2, function(x) all(x == 0)))
zeroes_X <- which(apply(X_FABIA, 2, function(x) all(x == 0)))
zeroes <- union(zeroes_B, zeroes_X)

if (length(zeroes) > 0) {
  X_FABIA <- as.matrix(X_FABIA[, -zeroes])
  B_FABIA <- as.matrix(B_FABIA[, -zeroes])
}


# PLOT
K_FABIA <- ncol(X_FABIA)

X_FABIA[X_FABIA > 10] <- 10
X_FABIA[X_FABIA < -10] <- -10


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

reorder <- union(c(13, 14,45, 16, 87, 84, 88, 89, 90, 91, 93, 96, 26,
                   23, 24,
                   20, 21, 22, 28, 47, 46, 51, 53, 56, 59, 60, 65, 69,94, 97,
                   74,73, 76, 77, 78, 79, 83,80, 85,
                   35, 59, 61, 69, 81,
                   36, 70, 47, 
                   39, 18, 33, 44,31,32,68,
                   43,  75, 82, 55, 37, 15, 27,38, 62, 98,40, 54,52,
                   58,42, 1, 4, 7, 11, 49, 95, 
                   34, 41, 48, 
                   2, 5,17, 25, 30, 63, 64, 8, 9,
                   19 ),
                 1:K_FABIA)

FABIA_plot <- plot_matrix(X_FABIA[ord,reorder ], title = "FABIA Factor Matrix")
FABIA_plot

pdf(file = "figures/zeisel_FABIA.pdf", width = 10, height = 8, pointsize = 16)

grid.arrange(group_plot, subclass_plot, FABIA_plot, nrow = 1, widths = c(5, 1, 25))

dev.off()


