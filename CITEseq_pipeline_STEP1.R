#PsychEncode fastq: /sc/arion/projects/CommonMind/roussp01a/Psychencode/MSSM_Aim3/fastq_files
#cellranger-count: /sc/arion/projects/CommonMind/wen/PEC_MSSM_snRNA/results
library(DropletUtils)
library(Seurat)
library(lattice)
library(XLConnect)
##Create directories
#source file from cDNA_dir and HTO_dir
######parameters to modify accordingly
cDNA_dir = "/sc/arion/projects/CommonMind/wen/PEC_MSSM_snRNA/results/"
HTO_dir = "/sc/arion/projects/CommonMind/hclee/shared/PEC_HTO_kite/"
NYGC_delivery = "ROU_14159_NAN_Lane_1"
#make sure to put NYGC manifest under work_dir
NYGC_manifest = "ROU_14159_NAN_Lane_1-LANESEQ-MTR-06370.xlsx"
work_dir = "/sc/arion/projects/CommonMind/shan/CITEseq_pipeline/test/"
DRY_RUN = T
QUEUE_PRIORITY = "premium"
ACCOUNT_FOR_JOBS = "acc_CommonMind"
QUEUE_STRING = NULL
PROJECT = "PEC_MSSM"
################
if (!file.exists(work_dir)) dir.create(work_dir, recursive = T)
script_dir = paste0(work_dir,"scripts/")
if (!file.exists(script_dir)) dir.create(script_dir, recursive = T)
demultiplexed_dir = paste0(work_dir,"HTO_demultiplexed/") 
if (!file.exists(demultiplexed_dir)) dir.create(demultiplexed_dir, recursive = T)
barcodeRanks_dir = paste0(work_dir,"barcodeRanks/")
if (!file.exists(barcodeRanks_dir)) dir.create(barcodeRanks_dir, recursive = T)
limitTsigF_dir = paste0(work_dir,"limitTsigF/")
if (!file.exists(limitTsigF_dir)) dir.create(limitTsigF_dir, recursive = T)
cellSNP_dir = paste0(work_dir,"cellSNP/")
if (!file.exists(cellSNP_dir)) dir.create(cellSNP_dir, recursive = T)
cellSNP_barcode_dir = paste0(cellSNP_dir,"barcodes/")
if (!file.exists(cellSNP_barcode_dir)) dir.create(cellSNP_barcode_dir, recursive = T)
log_dir = paste0(work_dir,"log/")
if (!file.exists(log_dir)) dir.create(log_dir, recursive = T)
###helper functions
submitJob = function(submitString,verbose=T){
    submitInfo=environment()
    submitInfo$submitString=submitString
    if(!DRY_RUN)
      submitInfo$systemReturnString=system(submitString, intern=T)
    else
      submitInfo$systemReturnString="Job <00000000> is submitted to queue <SOMETHING>"

    submitInfo$jobId=gsub("^.+<","",gsub(">.+$","",submitInfo$systemReturnString))
  
    if(verbose){
      print(submitInfo$submitString)
      print(submitInfo$systemReturnString)
      print(submitInfo$jobId)
    }
  
    return(submitInfo)
  }
  
write_csv = function(vec, fn){write.table(vec, file = fn, append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")}
##Read in all samples in a NYGC delivery seperated by pool ID
NYGC_wb = loadWorkbook(paste0(work_dir, NYGC_manifest),create = TRUE)
pool_data = readWorksheet(NYGC_wb, sheet = "Pool Contents", startRow = 6, endRow = 0, startCol = 2, endCol = 0, header = FALSE)
pool_data = pool_data[,1:2]
colnames(pool_data) = c("pool_ID","sample_ID")
#customize the screened strings
pool_data = pool_data[grepl("PEC_MSSM", pool_data[["sample_ID"]]) & grepl("cDNA", pool_data[["sample_ID"]]),] 
pool_data$sample_ID = gsub("_cDNA","",pool_data$sample_ID)
pool_data$sample_ID = gsub("_","-",pool_data$sample_ID)
pool_data$sample_ID = paste0("Sample_",pool_data$sample_ID)
pool_data_list = split(pool_data$sample_ID, f = pool_data$pool_ID)
##Create limitTsigF_df
limitTsigF_filename = paste0(limitTsigF_dir,NYGC_delivery,"_limitTsigF.csv")
write(file=limitTsigF_filename, paste0("sample_ID,limitTsigF_num"))

###The following loop execute DropletUtils swappedDrops(), EmptyDrops(), MULTIseqDemux() and HTODemux()
#####################################
for (j in 1:length(pool_data_list)) {
print(paste("process pool ID ", names(pool_data_list)[j]))
##pooled samples (read in all samples from NYGC manifest file in pipeline then divide)
#pooled_samples = c("Sample_PEC-MSSM-Set1","Sample_PEC-MSSM-Set4","Sample_PEC-MSSM-Set7")
pooled_samples = pool_data_list[[j]]
#read in cellranger aligned raw matrix
#sce = read10xCounts("/sc/arion/projects/CommonMind/wen/PEC_MSSM_snRNA/results/Sample_PEC-MSSM-Set1-cDNA/outs/raw_feature_bc_matrix")
#sce = Read10X("/sc/arion/projects/CommonMind/wen/PEC_MSSM_snRNA/results/Sample_PEC-MSSM-Set1-cDNA/outs/raw_feature_bc_matrix")
#Read in raw matrix and Remove swapped barcodes
set.seed(runif(1, min=0, max=10000))
mult.mol.info = paste0(cDNA_dir, pooled_samples, "-cDNA/outs/molecule_info.h5")
#mult.mol.info = c(paste0(cDNA_dir,"Sample_PEC-MSSM-Set1","-cDNA/outs/molecule_info.h5"),"/sc/arion/projects/CommonMind/wen/PEC_MSSM_snRNA/results/Sample_PEC-MSSM-Set4-cDNA/outs/molecule_info.h5","/sc/arion/projects/CommonMind/wen/PEC_MSSM_snRNA/results/Sample_PEC-MSSM-Set7-cDNA/outs/molecule_info.h5")
s.out <- swappedDrops(mult.mol.info, min.frac=0.9)
#limitTsigF = vector(mode = "integer", length = length(pooled_samples))

##Loop thru each sample in the pool
for (i in 1:length(pooled_samples)) {
print(paste("process ",pooled_samples[i]))
sce = s.out$cleaned[[i]]
#barcode ranks and plot
br.out = barcodeRanks(sce)
pdf(paste0(barcodeRanks_dir, pooled_samples[i] ,"_br_out.pdf"),width = 7, height = 7)
# Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))
dev.off()
#Remove empty drops
set.seed(runif(1, min=0, max=10000))
e.out <- emptyDrops(sce)
is.cell <- e.out$FDR <= 0.01
#sum(is.cell, na.rm=TRUE)
#sig cell table
#limitTsigF[i] = table(Limited=e.out$Limited, Significant=is.cell)[2,1]
write(file=limitTsigF_filename, paste0(pooled_samples[i],",",table(Limited=e.out$Limited, Significant=is.cell)[2,1]), append = TRUE)
sce = sce[,colnames(sce) %in% rownames(e.out)[e.out$FDR <= 0.01]]
##output barcodes for cellSNP
barcodes = data.frame(paste0(colnames(sce),"-1"))
colnames(barcodes) = "barcodes"
write.table(barcodes$barcodes, file = paste0(cellSNP_barcode_dir,pooled_samples[i],"_barcodes.tsv"), col.names = F, row.names = F, quote = F)
## find common cells in HTO and cDNA
## then run multiseq/htodemux to threshold hto counts
#gdo_count = Read10X(paste0(GDO_dir, '/counts/'), gene.column=1)
hto_count = Read10X(paste0(HTO_dir, pooled_samples[i] ,'/counts/'), gene.column=1)
#cells = read.csv(paste0(cDNA_dir, '/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'), header=FALSE)
#cells = substr(cells[,], 1, 16)
#hto_count = hto_count[ ,colnames(hto_count) %in% cells]
hto_count = hto_count[ ,colnames(hto_count) %in% colnames(sce)]

hto_obj <- CreateSeuratObject(counts = hto_count, assay = "HTO")
hto_obj <- NormalizeData(hto_obj, assay = "HTO", normalization.method = "CLR")
hto_obj = MULTIseqDemux(hto_obj, assay = "HTO", autoThresh=TRUE)
hto_obj = HTODemux(hto_obj, assay = "HTO")
    
write_csv(hto_obj$MULTI_ID, paste0(demultiplexed_dir,pooled_samples[i], '_MULTIseq.csv'))
write_csv(hto_obj$hash.ID, paste0(demultiplexed_dir,pooled_samples[i], '_HTOdemux.csv')) 
}
}
###############################
#cellSNP
cellSNP_script = paste0(script_dir,"cellSNP.sh")
write(file = cellSNP_script, paste0("#!/bin/bash
	#BSUB -J ", PROJECT, "
    #BSUB -q ", QUEUE_PRIORITY, "
    #BSUB -P ", ACCOUNT_FOR_JOBS, " 
    #BSUB -n 8 
    #BSUB -R 'span[hosts=1]'
    #BSUB -R 'rusage[mem=50000]'",
    QUEUE_STRING, "
    #BSUB -W 144:00 
    #BSUB -L /bin/bash
    #BSUB -oo ", log_dir, "cellSNP.%I.out
    #BSUB -eo ", log_dir, "cellSNP.%I.err
  
    module purge
    ml python
	
###############################
sample_id=$(sed \"${LSB_JOBINDEX}q;d\" <<< '", paste(pool_data$sample_ID, collapse="\n"), "')
cellSNP -s ",cDNA_dir, "${sample_id}-cDNA/outs/possorted_genome_bam.bam -b ", cellSNP_barcode_dir, "${sample_id}_barcodes.tsv -O ", cellSNP_dir, "${sample_id}_cellSNP_out -R /sc/arion/projects/CommonMind/shan/cellSNP_geno_ref/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz -p 20 --minMAF 0.1 --minCOUNT 20
"))

cellSNP_BatchJob = submitJob(   paste0("cat ", cellSNP_script, "| bsub -J '", PROJECT, "[1-",length(pool_data$sample_ID) ,"]'")   )

