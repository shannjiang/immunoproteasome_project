# Loading the package for proteasome analysis
library(DEP)
# Loading a package required for data handling
library(dplyr)
library(pheatmap)
library(ggplot2)
library(genefu)
library(plotly)
library(tidyr)

#helper function
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

work_dir = 'C:/Users/shann/Documents/Natura/Mayank_immunopeptidomics/'
prot_dat = read.csv(file = paste0(work_dir,'MK1479_Fusion_MSfragger_combined_peptide.csv'), header = T)
prot_dat[is.na(prot_dat)] = 0

# Are there any duplicated gene names?
prot_dat$Peptide.Sequence %>% duplicated() %>% any()

# Make a table of duplicated gene names
prot_dat %>% group_by(Gene) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

# Make unique names using the annotation in the "Gene" column as primary names and the annotation in "Peptide.Sequence" as name for those that do not have an gene name.
prot_unique <- make_unique(prot_dat,"Peptide.Sequence","Gene", delim = ";")

prot_sample_columns = c(63:86)
experimental_design = data.frame(label = c(paste0('LHA_',c(1:2)),'hTau_App_1','WT_1',paste0('LMP_',c(1:4)),paste0('App_',c(1:3)),'hTau_App_2','App_4',paste0('LMP_App_',c(1:4)),'WT_2','LHA_3','hTau_App_3','WT_3','LHA_4','hTau_App_4','WT_4'))
experimental_design$label2 = experimental_design$label
experimental_design = experimental_design %>% separate(label2,into=c('condition','replicate'),sep = -2,convert = TRUE)
experimental_design$replicate = gsub('\\_','',experimental_design$replicate)
colnames(prot_unique)[prot_sample_columns] = experimental_design$label
#bad sample removal
badSamples = c('App_2','hTau_App_1','hTau_App_2','LHA_1','LHA_2','WT_1','WT_2','WT_4')
prot_unique = prot_unique[,!colnames(prot_unique) %in% badSamples]
experimental_design = experimental_design[!experimental_design$label %in% badSamples,]
prot_sample_columns = c(63:78)

data_se <- make_se(prot_unique, prot_sample_columns, experimental_design)
data_se_parsed <- make_se_parse(prot_unique, prot_sample_columns)

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

##############filter out proteins of low quality
# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

################imputation
# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- DEP::impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- DEP::impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- DEP::impute(data_norm, fun = "knn", rowmax = 0.9)

# Plot intensity distributions before and after imputation data_imp
plot_imputation(data_norm, data_imp)

# Plot intensity distributions before and after imputation data_imp_man
plot_imputation(data_norm, data_imp_man)

# Plot intensity distributions before and after imputation data_imp_knn
plot_imputation(data_norm, data_imp_knn)

# Differential enrichment analysis  based on linear models and empherical Bayes statistics
# Test every sample versus control
data_diff <- test_diff(data_imp_knn, type = "all")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

#################visualization of the results
# Plot the first and second principal components (PCA)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 2, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

# Generate a results table
res <- get_results(dep)
write.csv(res, file = paste0(work_dir,'Mayank_immunopeptidomics_DE_res.csv'))

#pairwise heatmap
datExpr = as.data.frame(data_imp@assays@data)
datExpr = datExpr[,!colnames(datExpr) %in% c('group','group_name')]
#col_annot_df = data.frame(condition = c(rep('LHA',4),rep('hTau_APP',4),rep('LMP_APP',3),rep('APP',3),rep('LMP',3),rep('WT',3)))
col_annot_df = data.frame(condition = experimental_design$condition)
#col_annot_df$condition = gsub('batch2_','',experimental_design$label)
rownames(col_annot_df) = colnames(datExpr)
col_annot_df$sample = rownames(col_annot_df)

#sigGene heatmap helper function
pairwise_sigGene_heatmap = function(datExpr = datExpr, col_annot_df = col_annot_df, comparison = c('g1','g2'), res = res) {
col_annot_df2 = col_annot_df[with(col_annot_df,order(condition,sample)),]
col_annot_df2 = col_annot_df2[col_annot_df2$condition %in% comparison,]
col_annot_df2 = rbind(col_annot_df2[col_annot_df2$condition %in% comparison[1],],col_annot_df2[col_annot_df2$condition %in% comparison[2],])
datExpr2 = datExpr[,colnames(datExpr) %in% rownames(col_annot_df2)]
datExpr2 = datExpr2[,rownames(col_annot_df2)]
res4screen = res[,colnames(res) %in% c('name', paste0(comparison[1],'_vs_',comparison[2],'_p.val'), paste0(comparison[2],'_vs_',comparison[1],'_p.val'))]
sigRes = res4screen[res4screen[,2] < 0.05,]
datExpr2 = datExpr2[rownames(datExpr2) %in% sigRes$name,]
fig = pheatmap(datExpr2, scale = 'row', cluster_rows = T, cluster_cols = F, show_rownames = F, show_colnames = F, fontsize_row = 6, annotation_col = col_annot_df2, annotation_row = NA)
save_pheatmap_pdf(fig, paste0(work_dir,'Mayank_immunopeptidomics_',comparison[1],'_vs_',comparison[2],'_sigGene_heatmap.pdf'), width = 7, height = 7)
}

pairwise_sigGene_heatmap(datExpr,col_annot_df,c('App','hTau_App'),res)
pairwise_sigGene_heatmap(datExpr,col_annot_df,c('App','LMP_App'),res)
pairwise_sigGene_heatmap(datExpr,col_annot_df,c('hTau_App','LHA'),res)

#volcano plot helper function
#label_thre is the percentage of genes passing nominal significant 0.05
res4volcano = function(res = res, comparison = c('g1','g2'), label_thre = 0.75) {
	toMatch = c(paste0('^',comparison[1],'_vs_',comparison[2]),paste0('^',comparison[2],'_vs_',comparison[1]))
	colRetain = grep(paste(toMatch,collapse='|'),colnames(res))
	colRetain = c(c(1:2),colRetain)
	res2 = res[,colRetain]
	#get p.val colname
	pvalColNum = grep('p.val',colnames(res2))
	pvalColName = colnames(res2)[pvalColNum]
	res2['log10P'] = -log10(res2[pvalColName])
	res2['volcano_color'] = ifelse(res2[,pvalColName] < 0.05, 'black', 'grey')
	#label 75% of genes passed nominal significant 0.05
	nomSigRes = res2[res2$volcano_color == 'black',]
	threLog10P = quantile(nomSigRes$log10P, label_thre)
	res2$annot_thre = ifelse(-log10(res2[,pvalColName]) > threLog10P, 'TRUE', 'FALSE')
	return(res2)
}

#App_vs_hTau_App_ratio
res2 = res4volcano(res,c('App','hTau_App'),label_thre = 0.9)
groupVS = colnames(res2)[3]
groupVS = sub('_p.val','',groupVS)
groups = unlist(strsplit(groupVS,'_vs_'))
volcano_colr = as.character(unique(res2[,'volcano_color']))
res2$hTau_App_vs_App_ratio = -res2$App_vs_hTau_App_ratio
p = ggplot(res2,aes(x=hTau_App_vs_App_ratio,y=log10P)) + geom_point(aes(color=volcano_color, text = paste0('peptide: ',name,'\ngene: ',ID,'\nPadj: ',App_vs_hTau_App_p.adj,'\nratio: ',hTau_App_vs_App_ratio)),size=1.5,alpha=0.5) + scale_color_manual(breaks = unique(res2$volcano_color), values = volcano_colr) + theme(legend.position = "none") + xlab('Log2(fold change)') + ylab('-Log10(P value)') + geom_hline(yintercept=-log10(0.05),linetype='dashed',color='black') + annotate("text", x=min(res2[,6])+0.4, y=0, label= groups[1], size = 3) + annotate("text", x=max(res2[,6])-0.4, y=0, label= groups[2], size = 3)
ggplotly(p, tooltip = "text")

#App_vs_LMP_App_ratio
res2 = res4volcano(res,c('App','LMP_App'),label_thre = 0.9)
groupVS = colnames(res2)[3]
groupVS = sub('_p.val','',groupVS)
groups = unlist(strsplit(groupVS,'_vs_'))
volcano_colr = as.character(unique(res2[,'volcano_color']))
res2$LMP_App_vs_App_ratio = -res2$App_vs_LMP_App_ratio
p = ggplot(res2,aes(x=LMP_App_vs_App_ratio,y=log10P)) + geom_point(aes(color=volcano_color, text = paste0('peptide: ',name,'\ngene: ',ID,'\nPadj: ',App_vs_LMP_App_p.adj,'\nratio: ',LMP_App_vs_App_ratio)),size=1.5,alpha=0.5) + scale_color_manual(breaks = unique(res2$volcano_color), values = volcano_colr) + theme(legend.position = "none") + xlab('Log2(fold change)') + ylab('-Log10(P value)') + geom_hline(yintercept=-log10(0.05),linetype='dashed',color='black') + annotate("text", x=min(res2[,6])+0.4, y=0, label= groups[1], size = 3) + annotate("text", x=max(res2[,6])-0.4, y=0, label= groups[2], size = 3)
ggplotly(p, tooltip = "text")

#LHA_vs_hTau_App_ratio
res2 = res4volcano(res,c('hTau_App','LHA'),label_thre = 0.9)
groupVS = colnames(res2)[3]
groupVS = sub('_p.val','',groupVS)
groups = unlist(strsplit(groupVS,'_vs_'))
volcano_colr = as.character(unique(res2[,'volcano_color']))
#res2$hTau_App_vs_LHA_ratio = -res2$LHA_vs_hTau_App_ratio
p = ggplot(res2,aes(x=LHA_vs_hTau_App_ratio,y=log10P)) + geom_point(aes(color=volcano_color, text = paste0('peptide: ',name,'\ngene: ',ID,'\nPadj: ',LHA_vs_hTau_App_p.adj,'\nratio: ',LHA_vs_hTau_App_ratio)),size=1.5,alpha=0.5) + scale_color_manual(breaks = unique(res2$volcano_color), values = volcano_colr) + theme(legend.position = "none") + xlab('Log2(fold change)') + ylab('-Log10(P value)') + geom_hline(yintercept=-log10(0.05),linetype='dashed',color='black') + annotate("text", x=min(res2[,6])+0.4, y=0, label= groups[2], size = 3) + annotate("text", x=max(res2[,6])-0.4, y=0, label= groups[1], size = 3)
ggplotly(p, tooltip = "text")



#plotly
fig <- plot_ly(type = 'scatter', mode = 'markers') 
fig <- fig %>%
  add_trace(
    x = c(1:5), 
    y = rnorm(5, mean = 5),
    text = c("Text A", "Text B", "Text C", "Text D", "Text E"),
    hoverinfo = 'text',
    marker = list(color='green'),
    showlegend = F
  )

#ggplotly
geom_jitter(aes(text = paste0("Name: ", name, "\nYear: ", year2, "\nCount of issues: ", Scope)))
ggplotly(p, tooltip = "text")