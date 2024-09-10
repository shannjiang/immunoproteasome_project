library(dplyr)
library(ggplot2)
library(WebGestaltR)
library(cowplot)
library(DESeq2)

###dejager rosmap
dejager_pseudobulk_dir = '/home/shann/Documents/synapse/ROSMAP/snrnaseq/DeJager/pseudoBulk/'
out_dir = '/home/shann/Documents/Natura/immunoproteasome_project/human_snRNAseq/dejager_rosmap/IP_linear_regression/'
if(!dir.exists(out_dir)){dir.create(out_dir,recursive = T)}

proteasome_gene_id_df = read.table(file = paste0('/home/shann/Documents/Natura/proteosome/gene_id.txt'), header = T, sep = '\t')
#immunoproteasome_genes: PSMB8 PSMB9 PSMB10
proteasome_gene_id_df = proteasome_gene_id_df[!proteasome_gene_id_df$OGS %in% c('NFE2L1','NFE2L2','ADRM1','MAFF','MAFG','MAFK'),]
IP_gene_id_df = proteasome_gene_id_df[proteasome_gene_id_df$OGS %in% c('PSMB8','PSMB9','PSMB10'),]

#parameters for GSEA
enrichFullNames = c('geneontology_Biological_Process_noRedundant','geneontology_Cellular_Component_noRedundant','geneontology_Molecular_Function_noRedundant','pathway_KEGG')
enrichShortNames = c('BP','CC','MF','KEGG')

#helper
IP_Expr2meta = function(vst = msbb_bm10_vst, meta = msbb_bm10_meta){
  ip_vst = vst[rownames(vst) %in% IP_gene_id_df$OGS,]
  ip_vst = as.data.frame(t(ip_vst))
  meta$ip_expr = rowMeans(ip_vst)
  return(meta)
}

linear_regression = function(vst = msbb_bm10_vst, meta = msbb_bm10_meta, variates = variates, out_dir = out_dir, filename = filename){
  rownames(vst) = gsub('\\-','_',rownames(vst))
  rownames(vst) = gsub('\\/','__',rownames(vst))
  vst_meta = cbind(meta,as.data.frame(t(vst)))
  genes = rownames(vst)
  Estimate = vector()
  SE = vector()
  t_value = vector()
  P_value = vector()
  for(gene in genes) {
    formula = paste0(gene,'~',paste(variates,collapse = '+'))
    m = lm(formula,data = vst_meta)
    Estimate[length(Estimate) + 1] = coef(summary(m))[2,1]
    SE[length(SE) + 1] = coef(summary(m))[2,2]
    t_value[length(t_value) + 1] = coef(summary(m))[2,3]
    P_value[length(P_value) + 1] = coef(summary(m))[2,4]
  }
  genes = gsub('__','/',genes)
  genes = gsub('_','-',genes)
  res = data.frame(gene = genes, Estimate = Estimate, SE = SE, t_value = t_value, P_value = P_value)
  res$BH = p.adjust(res$P_value, method='BH', n = length(genes))
  res = res[order(res$P_value),]
  write.csv(res,file=paste0(out_dir,filename))
}

celltypes = c('excitatory','inhibitory','microglia','astrocytes','oligodendrocytes','opcs','endo')
covariates = c('ip_expr','sex','age_death','pmi','sizeFactor')

for(celltype in celltypes) {
  paste0('process ',celltype)
  z = load(file = paste0(dejager_pseudobulk_dir,celltype,'_pseudoBulk.RData'))
  #DESeq2 vst
  pseudoBulk_cts = as.data.frame(t(pseudoBulk_cts))
  pseudoBulk_meta$sex = as.factor(pseudoBulk_meta$sex)
  pseudoBulk_meta$pAD = as.factor(pseudoBulk_meta$pAD)
  #pseudoBulk_meta$braaksc = as.integer(pseudoBulk_meta$braaksc)
  expr_dds = DESeqDataSetFromMatrix(countData = pseudoBulk_cts, colData = pseudoBulk_meta,design = ~ -1)
  design(expr_dds) = formula(~ sex + age_death + pmi + pAD)
  pseudoBulk_vst = vst(expr_dds, blind = T)
  pseudoBulk_vst_assay = as.data.frame(assay(pseudoBulk_vst))
  meta = as.data.frame(colData(pseudoBulk_vst))
  pseudoBulk_meta = IP_Expr2meta(pseudoBulk_vst_assay,meta)
  pseudoBulk_vst = pseudoBulk_vst_assay
  linear_regression(pseudoBulk_vst,pseudoBulk_meta,covariates,out_dir,paste0('dejager_rosmap_snrnaseq_',celltype,'_IPexpr_on_gene_adj4SexAgedeathPMIsizeFactor.csv'))
}

#rnk GSEA
#turn off scientific notation
options(scipen=999)
gsea_dir = paste0(out_dir,'GSEA/')
if(!dir.exists(gsea_dir)){dir.create(gsea_dir,recursive=T)}

#generate rnk file
#res = read.csv(file = paste0(out_dir,'msbb_bm10_IPexpr_on_gene_adj4RINSexRaceAgedeathPMIsizeFactor.csv'), header = T, row.names = 1)
#res$rnk = -res$t_value
#rnk_df = res[,colnames(res) %in% c('gene','rnk')]
#write.table(rnk_df,file = paste0(out_dir,'gsea.rnk'),quote = F, row.names = F, col.names = F, sep = '\t')

#helper
rnk_GSEA = function(rnk = paste0(out_dir,'gsea.rnk'), enrichFullNames = enrichFullNames, enrichShortNames = enrichShortNames, geneIDtype = 'genesymbol', name = 'msbb_bm10', dir = gsea_dir){
  outcome_enrichment_dir = paste0(dir,name,'/')
  if(!dir.exists(outcome_enrichment_dir)){dir.create(outcome_enrichment_dir,recursive=T)}
  for (i in 1:length(enrichFullNames)) {
    enrichResult <- try(WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
      enrichDatabase=enrichFullNames[i], interestGeneFile=rnk,
      interestGeneType=geneIDtype, referenceGeneFile=refFile,
      referenceGeneType=geneIDtype, isOutput=TRUE,
      outputDirectory=outcome_enrichment_dir, projectName=paste0(name,'_',enrichShortNames[i]), fdrThr=1))
  }
}

#create a dataset vector and res file vector
#datasets = c('msbb_bm22','msbb_bm36','msbb_bm44','rosmap_dlpfc','mayo_tcx','mayo_cer')
#resFiles = c('msbb_bm22_IPexpr_on_gene_adj4RINSexRaceAgedeathPMIsizeFactor.csv','msbb_bm36_IPexpr_on_gene_adj4RINSexRaceAgedeathPMIsizeFactor.csv','msbb_bm44_IPexpr_on_gene_adj4RINSexRaceAgedeathPMIsizeFactor.csv','rosmap_dlpfc_IPexpr_on_gene_adj4sexEducRaceAgedeathPMIRINsizeFactor.csv','mayo_tcx_IPexpr_on_gene_adj4sexAgedeathRINsizeFactor.csv','mayo_cer_IPexpr_on_gene_adj4sexAgedeathRINsizeFactor.csv')

datasets = c('excitatory','inhibitory','microglia','astrocytes','oligodendrocytes','opcs','endo')
resFiles = c('dejager_rosmap_snrnaseq_excitatory_IPexpr_on_gene_adj4SexAgedeathPMIsizeFactor.csv','dejager_rosmap_snrnaseq_inhibitory_IPexpr_on_gene_adj4SexAgedeathPMIsizeFactor.csv','dejager_rosmap_snrnaseq_microglia_IPexpr_on_gene_adj4SexAgedeathPMIsizeFactor.csv','dejager_rosmap_snrnaseq_astrocytes_IPexpr_on_gene_adj4SexAgedeathPMIsizeFactor.csv','dejager_rosmap_snrnaseq_oligodendrocytes_IPexpr_on_gene_adj4SexAgedeathPMIsizeFactor.csv','dejager_rosmap_snrnaseq_opcs_IPexpr_on_gene_adj4SexAgedeathPMIsizeFactor.csv','dejager_rosmap_snrnaseq_endo_IPexpr_on_gene_adj4SexAgedeathPMIsizeFactor.csv')

for(i in 1:length(datasets)){
  print(paste0('process ',datasets[i]))
  res = read.csv(file = paste0(out_dir,resFiles[i]), header = T, row.names = 1)
  res$rnk = -res$t_value
  rnk_df = res[,colnames(res) %in% c('gene','rnk')]
  write.table(rnk_df,file = paste0(out_dir,'gsea.rnk'),quote = F, row.names = F, col.names = F, sep = '\t')
  rnk_GSEA(rnk = paste0(out_dir,'gsea.rnk'), enrichFullNames = enrichFullNames, enrichShortNames = enrichShortNames, geneIDtype = 'genesymbol', name = datasets[i], dir = gsea_dir)
}

#helper functions for dot plot by dataset
maxPower = function(vector){
    vector = vector[vector != 0]
    return(ceiling(-log10(vector[1])))
}

enrichmentTopDf = function(projectName = projectName, enrichmentDir = enrichmentDir, top = 5){
bp_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_BP/enrichment_results_',projectName,'_BP.txt'), header = T, sep = '\t')
cc_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_CC/enrichment_results_',projectName,'_CC.txt'), header = T, sep = '\t')
mf_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_MF/enrichment_results_',projectName,'_MF.txt'), header = T, sep = '\t')
kegg_df = read.table(file = paste0(enrichmentDir,'Project_',projectName,'_KEGG/enrichment_results_',projectName,'_KEGG.txt'), header = T, sep = '\t')

bp_df$FDR = gsub('^0$',paste0('1e-',maxPower(bp_df$FDR)),bp_df$FDR)
bp_df$FDR = as.numeric(bp_df$FDR)
cc_df$FDR = gsub('^0$',paste0('1e-',maxPower(cc_df$FDR)),cc_df$FDR)
cc_df$FDR = as.numeric(cc_df$FDR)
mf_df$FDR = gsub('^0$',paste0('1e-',maxPower(mf_df$FDR)),mf_df$FDR)
mf_df$FDR = as.numeric(mf_df$FDR)
kegg_df$FDR = gsub('^0$',paste0('1e-',maxPower(kegg_df$FDR)),kegg_df$FDR)
kegg_df$FDR = as.numeric(kegg_df$FDR)
bp_df$log10fdr = -log10(bp_df$FDR)
cc_df$log10fdr = -log10(cc_df$FDR)
mf_df$log10fdr = -log10(mf_df$FDR)
kegg_df$log10fdr = -log10(kegg_df$FDR)

#add category
bp_df$cat = 'BP'
cc_df$cat = 'CC'
mf_df$cat = 'MF'
kegg_df$cat = 'KEGG'

#combine the top five terms from each category
top_df = rbind(bp_df[1:top,],cc_df[1:top,],mf_df[1:top,],kegg_df[1:top,])
top_df$description = factor(top_df$description, levels = rev(top_df$description))
return(top_df)
}

enrichmentDotPlot = function(top_df) {
p = ggplot(top_df, aes(x=enrichmentScore, y = description, color = log10fdr, size = size)) +
  geom_point() +
  scale_colour_gradient2(name = '-log10(FDR)', low = 'blue', mid = 'grey', high = 'red', midpoint = midpoint) +
  cowplot::theme_cowplot() +
  theme(axis.text = element_text(size = 8, angle = 0, vjust = 1, hjust=1), axis.title=element_text(size=8), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  ylab('') +
  theme(axis.ticks = element_blank()) + geom_hline(yintercept = c(5.5,10.5,15.5),linetype = 'dashed')
#dev.off()
return(p)
}

#datasets = c('msbb_bm10','msbb_bm22','msbb_bm36','msbb_bm44','rosmap_dlpfc','mayo_tcx','mayo_cer')
#dotplot by dataset
for (dataset in datasets) {
  enrichmentDir = paste0(gsea_dir,dataset,'/')
  #projectName = outcomes[1]
  print(paste0('process ',dataset))
  top_df = enrichmentTopDf(dataset,enrichmentDir,5)
  midpoint = mean(top_df$log10fdr)
  p = enrichmentDotPlot(top_df)
  print(p)
  pdf(file = paste0(enrichmentDir, dataset,'_top_enrichment_by_category.pdf'), width = 7, height = 0.2*length(top_df$description))
  print(p)
  dev.off()
}

#top GO Term enrichment across datasets
datasetTopGOTermEnrichmentGraph = function(datasets = datasets, out_dir = gsea_dir, top = 1) {
  enrichList = list()
  GOshortNames = c('BP','CC','MF')
  for (dataset in datasets) {
    for (GOcat in GOshortNames) {
      if (file.exists(paste0(out_dir,dataset,'/Project_',dataset,'_',GOcat,'/enrichment_results_',dataset,'_',GOcat,'.txt'))) {
        df = read.table(file = paste0(out_dir,dataset,'/Project_',dataset,'_',GOcat,'/enrichment_results_',dataset,'_',GOcat,'.txt'), sep = '\t', header = T)
        df$dataset = dataset
        enrichList[[length(enrichList)+1]] = df
      }
    }
  }
  #keep the top term(s) in each GO category and dataset
  braakEnrichTop = do.call(rbind, (lapply(enrichList, function(x) x[top,])))
  braakEnrichTopAll = do.call(rbind, (lapply(enrichList, function(x){x[x$description %in% braakEnrichTop$description,]})))
  #remove description with too long name
  braakEnrichTopAll$description = gsub(',.*','',braakEnrichTopAll$description)
  braakEnrichTopAll$FDR = gsub('^0$',paste0('1e-',maxPower(braakEnrichTopAll$FDR)),braakEnrichTopAll$FDR)
  braakEnrichTopAll$FDR = as.numeric(braakEnrichTopAll$FDR)
  braakEnrichTopAll$log10p = -log10(braakEnrichTopAll$FDR)

  braakEnrichTopAll$dataset = factor(braakEnrichTopAll$dataset,levels = datasets)
  braakEnrichTopAll4barplot = braakEnrichTopAll[c('description','normalizedEnrichmentScore','dataset','log10p')]
  ###arrange the order of description with mitochondrial related GO terms on the top
  #mitoTerms = c('NADH dehydrogenase complex assembly','NADH dehydrogenase complex','mitochondrial protein complex','mitochondrial membrane part')
  #mitoTermsIn = mitoTerms[mitoTerms %in% braakEnrichTopAll4barplot$description]
  #braakEnrichTopNoMito = braakEnrichTopAll4barplot[!braakEnrichTopAll4barplot$description %in% mitoTermsIn,]
  #braakEnrichTopNoMito = braakEnrichTopNoMito[,colnames(braakEnrichTopNoMito) %in% c('description','normalizedEnrichmentScore')]
  #braakEnrichTopNoMitoAgg = aggregate(.~description,braakEnrichTopNoMito,mean)
  #braakEnrichTopNoMitoAgg = braakEnrichTopNoMitoAgg[order(braakEnrichTopNoMitoAgg$normalizedEnrichmentScore,decreasing = TRUE),]
  #descriptionLevels = c(braakEnrichTopNoMitoAgg$description,mitoTermsIn)
  #braakEnrichTopAll4barplot$description = factor(braakEnrichTopAll4barplot$description,levels = descriptionLevels)
  #pdf(file = paste0(out_dir,outcome,'_topGOterm.pdf'), width = 35, height = 5)
  p = ggplot(data = braakEnrichTopAll4barplot, aes(x = log10p, y = description, fill = normalizedEnrichmentScore)) + geom_bar(stat = 'identity') + scale_fill_gradient(low="darkblue", high="red") + geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'black') + facet_wrap(dataset~., strip.position='top',nrow = 1) + xlab('-log10(p-value)') + ylab('') + theme(text = element_text(size=20), axis.text.y = element_text(angle = 0),legend.position="bottom")
  #dev.off()
  return(p)
}

top = 3
p = datasetTopGOTermEnrichmentGraph(datasets = datasets, out_dir = gsea_dir, top = top)
print(p)
pdf(file = paste0(gsea_dir,'top',top,'_GOterm.pdf'), width = 35, height = 5*top)
print(p)
dev.off()
