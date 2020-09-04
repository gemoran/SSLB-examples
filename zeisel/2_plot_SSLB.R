# R script to get SSLB plots for Zeisel et al (2015) data + see whether marker genes are in biclusters

source('../SSLB_functions.R', echo=TRUE)

library(pals)

load("SSLB_result/SSLB_output.RData")

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
  


X_SSLB <- out$X
X_SSLB[X_SSLB > 10] <- 10
X_SSLB[X_SSLB < -10] <- -10

B_SSLB <- out$B
K_SSLB <- out$K


reorder <- union(c(4, 10,
                   5, 27, 26, 19, 29, 24,
                   3, 16, 21, 33,
                   1, 7, 14, 6, 22, 23,
                   2, 28, 15, 20, 11, 40, 
                   32, 13, 39, 42, 41,
                   8, 34, 31, 30, 25,
                   9, 18, 12,
                   35, 43, 55, 63,
                   17, 38, 43), 1:out$K)

X_SSLB <- X_SSLB[, reorder]
B_SSLB <- B_SSLB[, reorder]

X_plot <- plot_matrix(X_SSLB[ord, ], title = "SSLB Factor Matrix") 
  
X_plot


pdf(file = "figures/zeisel_SSLB_all.pdf", width = 10, height = 8, pointsize = 16)
grid.arrange(group_plot, subclass_plot, X_plot, nrow = 1, widths = c(5, 1, 25))
dev.off()


#--------------------------
# in group 2, SSLB has grouped (none) into another category, let's investigate

# subclasses in group 2
labels <- unique(meta_data$level2class[meta_data$group == 2])

test <- data.frame(group = groups, sc = subclasses, X = X_SSLB[,10])
test <- test[ord, ]

n_sub <- length(labels)
y_pos <- numeric(n_sub)
tot <- sum(groups == 2)
temp <- sum(meta_data$level2class[meta_data$group == 2] == labels[1] )
y_pos[1] <- tot - temp/2
old <- tot - temp

for (i in 2:n_sub) {
  temp <- sum(meta_data$level2class[meta_data$group == 2] == labels[i] )
  y_pos[i] <- old - temp/2
  old <- old - temp
}

label_data <- data.frame(x = rep(-1, n_sub), 
                         y = y_pos,
                         labels = labels)

cols <- brewer.paired(n_sub)
cols[1] <- "black"

g1 <- plot_matrix(test$sc[test$group==2], title = "") + 
  scale_fill_gradientn(colors = cols) +
  theme(plot.margin=unit(c(5.5, 5.5, 10, 30),"pt")) +
  geom_text(data = label_data, aes(x = x, y = y, label = labels), size = 4.5, angle = 160, vjust = 1) +
  coord_cartesian(clip = 'off') +
  theme(panel.border =element_blank()) 

g2 <- plot_matrix(test$X[test$group == 2], title = "") +
  theme(plot.margin=unit(c(5.5, 5.5, 10, 5.5),"pt")) 

pdf("figures/group2_zoom_in.pdf", width = 1.75, height = 8.5)
grid.arrange(g1, g2, nrow = 1, widths = c(4, 1.5))
dev.off()

# checking for gene markers identified by Zeisel
# interneurons
which(B_SSLB[which(genes == "Pnoc"), ] !=0)

# s1 pyramidal
which(B_SSLB[which(genes == "Tbr1"), ] !=0)
which(B_SSLB[which(genes == "Gm11549"), ] !=0)

# CA pyramidal
which(B_SSLB[which(genes == "Spink8"), ] !=0)

#oligodendrocytes
which(B_SSLB[which(genes == "Hapln2"), ] !=0)


#Hap1n2 is in bicluster 44 - what is there?
sum(X_SSLB[,44]!=0)
sum(X_SSLB[groups == 4,44]!=0)
# only oligodendrocytes! 
sc <- plot_matrix(subclasses[ord[groups == 4]]) + scale_fill_gradientn(colors = brewer.paired(length(unique(subclasses[groups == 4]))))
test <- X_SSLB[,44]
test[test!=0] <- 1
test_plot <- plot_matrix(test[ord[groups == 4]])
grid.arrange(sc, test_plot, nrow = 1)
# but not corresponding to any Zeisel oligodendrocyte subtype..........



# endothelial
which(B_SSLB[which(genes == "Ly6c1"), ] !=0)

# Ly6c1 is in bicluster 49 
sum(X_SSLB[,49]!=0)
sum(X_SSLB[groups == 6,49]!=0)
sum(X_SSLB[groups == 7,49]!=0)
# 7 endothelial, 1 astrocyte

# Ly6c1 also in bicluster 76
sum(X_SSLB[, 76]!=0)
sum(X_SSLB[groups == 6, 76]!=0)
sum(X_SSLB[groups == 5, 76]!=0)
# 4 endothelial,1 microglia

# mural
which(B_SSLB[which(genes == "Acta2"), ] !=0)


sum(X_SSLB[, 33]!=0)

sum(X_SSLB[groups == 4, 33]!=0)
sum(X_SSLB[groups == 5, 33]!=0)
sum(X_SSLB[groups == 6, 33]!=0)
sum(X_SSLB[groups == 7, 33]!=0)
sum(X_SSLB[groups == 9, 33]!=0)


sum(X_SSLB[, 81]!=0)
sum(X_SSLB[groups == 9, 33]!=0)
sum(X_SSLB[groups == 7, 33]!=0)

# How many subtypes did SSLB find?
# interneurons
length(unique(subclasses[groups == 1]))
apply(X_SSLB[groups == 1, ], 2, function(x) sum(x!=0))

# s1 pyramid
length(unique(subclasses[groups == 2]))

apply(X_SSLB[groups == 2, ], 2, function(x) sum(x!=0))

#------------------------------------------

# residual plot ----------------
Y = as.matrix(read.table("data/Y_quantile.txt", stringsAsFactors = F))

residuals = Y - X_SSLB %*% t(B_SSLB)
dat = melt(residuals)
dat$y = as.vector(Y)
dat$y_fit = as.vector(X_SSLB %*% t(B_SSLB))

dat_filter = dat[abs(dat$value) < 10, ]

pdf("figures/SSLB_residual_hist_filter.pdf", width = 4, height = 3.2)

ggplot(dat_filter, aes(x = value)) + 
  geom_histogram(bins = 50) +
  labs(title = "Zeisel Data: SSLB Residuals")  

dev.off()

pdf("figures/SSLB_residual_hist.pdf", width = 4, height = 3.2)

ggplot(dat, aes(x = value)) + 
  geom_histogram(bins = 50) +
  labs(title = "Zeisel Data: SSLB Residuals")  

dev.off()

dat_filter = dat
dat_filter = dat_filter[order(sample(1:nrow(dat_filter), 10000, replace = F)), ]
dat_filter$n = 1:10000

pdf("figures/SSLB_fitted_plot.pdf", width = 4, height = 3.2)

ggplot(dat_filter, aes(x = y, y = y_fit)) + 
  geom_point() +
  labs(title = "Zeisel Data: SSLB") + 
  ylab("Y fitted") +
  xlab("Y") + 
  geom_abline(intercept = 0, slope = 1, color = "red")

dev.off()

dat_filter = dat
dat_filter = dat_filter[order(sample(1:nrow(dat_filter), 10000, replace = F)), ]
dat_filter$n = 1:10000

dat_filter$y_fit[dat_filter$y_fit < 0] = 0

pdf("figures/SSLB_fitted_log_plot.pdf", width = 4, height = 3.2)

ggplot(dat_filter, aes(x = log(y + 1e-5), y = log(y_fit + 1e-5))) + 
  geom_point() +
  labs(title = "Zeisel Data: SSLB") + 
  ylab("log(Y fitted)") +
  xlab("log(Y)") + 
  geom_abline(intercept = 0, slope = 1, color = "red")

dev.off()

#-------
# range of sparsity
sparse_dat = data.frame(X = apply(X_SSLB, 2, function(x) sum(x!=0))/N * 100,
                        B = apply(B_SSLB, 2, function(x) sum(x!=0))/G * 100,
                        n = factor(1:K_SSLB))
sparse_dat = melt(sparse_dat)

pdf("figures/sparsity_levels.pdf", width = 13, height = 3)

ggplot(sparse_dat, aes(x = n, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_y_continuous(limits = c(0, 100)) + 
  ylab("% non-zero") + 
  xlab("Bicluster") +
  labs(title = "Zeisel Data: SSLB Sparsity Levels")

dev.off()
