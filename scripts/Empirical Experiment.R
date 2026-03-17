library(MASS)
library(copula)
library(Matrix)
library(ggplot2)

# load Huntley counts
X_real <- read.table("Huntley_counts.txt", header = TRUE, row.names = NULL)

# transpose if genes are rows
X_real <- t(as.matrix(X_real))

dim(X_real)   # samples × genes

# keep top variable genes
vars <- apply(X_real,2,var)

top_genes <- order(vars, decreasing=TRUE)[1:200]

X_real <- X_real[,top_genes]

#Log transform
X_real <- log1p(X_real)

###Fit Gaussian Copula model 

# estimate covariance
Sigma_real <- cov(X_real)

# sample synthetic data
n <- nrow(X_real)

X_syn <- mvrnorm(n, mu = colMeans(X_real), Sigma = Sigma_real)


#Diagnostic A: Marginal fidelity 

###KS test for each gene
ks_values <- sapply(1:ncol(X_real), function(j){
  ks.test(X_real[,j], X_syn[,j])$statistic
})

mean_ks <- mean(ks_values)
mean_ks


##Plot 
ks_df <- data.frame(KS = ks_values)

ggplot(ks_df, aes(KS)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "darkblue") +
  labs(
    title = "Distribution of KS distances across genes (Alzheimer's dataset)",
    x = "Kolmogorov–Smirnov distance",
    y = "Number of genes"
  ) +
  theme_minimal()

#Diagnostic B: Covariance fidelity 

Sigma_syn <- cov(X_syn)

cov_diff <- norm(Sigma_real - Sigma_syn, type="F")

cov_diff

##Heatmap
library(pheatmap)

pheatmap(Sigma_real - Sigma_syn,
         show_rownames=FALSE,
         show_colnames=FALSE)

pheatmap(
  Sigma_real - Sigma_syn,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Covariance difference: Real vs Synthetic (Alzheimer's dataset)",
  legend = TRUE
)

pheatmap(
  Sigma_real - Sigma_syn,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Covariance difference: Real vs Synthetic (Alzheimer's dataset)",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  legend = TRUE
)

#covariance difference matrix 
#Blue → negative difference
#	White → near 0 (good match)
#	Red → positive difference.



#Diagnostic C: Downstream regression

# create synthetic phenotype
y <- rowMeans(X_real[,1:10]) + rnorm(nrow(X_real),0,0.5)


##Fit regression model

fit_real <- lm(y ~ X_real)
fit_syn  <- lm(y ~ X_syn)

beta_real <- coef(fit_real)[-1]
beta_syn  <- coef(fit_syn)[-1]


###Compare coefficients

df_beta <- data.frame(beta_real, beta_syn)


ggplot(df_beta, aes(beta_real, beta_syn)) +
  geom_point() +
  geom_abline(color = "red") +
  labs(
    title = "Regression coefficient comparison: Real vs Synthetic (Alzheimer's dataset)",
    x = "Regression coefficients (Real data)",
    y = "Regression coefficients (Synthetic data)"
  ) +
  theme_minimal()