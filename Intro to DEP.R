source("https://bioconductor.org/biocLite.R")
biocLite("DEP")

library("DEP")
library("dplyr")

#For LFQ - Maz Quant analysis - Shiny web app
run_app("LFQ")

#Loading data
data <- UbiLength
data <- filter(data, Reverse != "+", Potential.contaminant != "+")
dim(data)
colnames(data)

#Duplicate genes names?
data$Gene.names %>% duplicated() %>% any()
#Make table displaying duplicated genes
data %>% group_by(Gene.names) %>% summarise(frequency = n()) %>% arrange(desc(frequency)) %>%
  filter(frequency > 1)

data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

data_unique$name %>% duplicated() %>% any()

#Generate a SummarizedExperiment object
LFQ_columns <- grep("LFQ.", colnames(data_unique))
experimental_design <- UbiLength_ExpDesign
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

data_se_parsed <- make_se_parse(data_unique, LFQ_columns)
data_se

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

#Filter out proteins containing too many missing data points
data_filt <- filter_missval(data_se, thr = 0)
data_filt2 <- filter_missval(data_se, thr=1)

plot_numbers(data_filt)
plot_coverage(data_filt)

#Normalization
data_norm <- normalize_vsn(data_filt)
meanSdPlot(data_norm)
#See how data is distributed following normalization
plot_normalization(data_filt, data_norm)
#Heatmap of proteins with missing values
#Use this to investigate the nature of the missing values
plot_missval(data_filt)

#plot distributions and cumulative fraction of proteins with and wihtout missing values
#Here we see that the missing values are enriched in low count proteins, 
#suggesting that they may be below the detection level in many samples. 
#Combined, the heatmap and the dist. plot suggest that the values are missing, not at random (MNAR)
plot_detect(data_filt)

data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

plot_imputation(data_norm, data_imp)
plot_imputation(data_norm, data_imp_man)

#Differential enrichment analysis
data_diff <- test_diff(data_imp, type = "control", control = "Ctrl")
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

#Define cutoffs to make significant proteins
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

#Visualize results

plot_pca(dep, x=1, y=2, n=500, point_size = 4)

#correlation matrix - correlation between samples
plot_cor(dep, significant = T, lower = 0, upper = 1, pal = "Reds")

#Heatmap of sig proteins
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))
#Volcano plot of specifi contrasts
plot_volcano(dep, contrast = "Ubi6_vs_Ctrl", label_size = 2, add_names = T)

#barplots for proteins of interest
plot_single(dep, proteins = c("USP15", "IKBKG"))
plot_single(dep, proteins = "USP15", type = "centered")

# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)

# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

colnames(data_results)

sessionInfo()
