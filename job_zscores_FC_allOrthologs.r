#!/Local/md_shenorr/R-3.2.0/bin/Rscript
#script is ment to be run on ATLAS one job per file. arguments recieved: RDSfile,disease,RFST
#goal- take .data files of microarray files, and calculate zscores and FC for each file WITH ALL ORTHOLOGS! INCLUDING ONE2MANY, MANY2MANY.

require(GEOdb)
require(applyBy)

# RFST="RF"
# disease="Influenza"
# RDSfile="HS_EpilepsyAutism_1_data.rds"

disease <- Sys.getenv("disease")
RDSfile <- Sys.getenv("RDSfile")
RFST <- Sys.getenv("RFST")


#map probe to gene
getMapping = function(eset)
{
  message(" *** In getMapping *** ")
  GPL = annotation(eset)
  annot <- paste(substr(alibrary(GPL), 1, nchar(alibrary(GPL))-3))
  mapping = eval(parse(text = paste0(annot, "ENTREZID")))
  return(mapping)
}

data_path = paste0("/storage/md_shenorr/rachelly/Shared_Dir/GS_analysis/Microarrays/",RFST,"/",disease,"/V1.9/")
if (disease=="Sepsis"){
  data_path = "/storage/md_shenorr/rachelly/Shared_Dir/GS_analysis/Microarrays/ST/Sepsis/V1.9_8best/" #todo validate with rachely!!
}
setwd(data_path)
eset = readRDS(RDSfile)
#if (is.null(probe_eset)) stop("Error: probe_eset is NULL, this is probably a gene-level eset - need to adjust code for this in CalcZscores")
mapping = getMapping(eset)                  
p2genes <- unlist(as.list(mapping[mappedkeys(mapping)])) # probe->gene mapping
entrez_features <- rownames(eset) %in% names(p2genes)    # finding Entrez IDs in mapping
annotated_probe_eset <- eset[entrez_features,]           # restricting the eset to mapped probes only
probe_factor <- p2genes[rownames(annotated_probe_eset)]
mean_per_probe = rowMeans(exprs(annotated_probe_eset),na.rm = T)   # calculating the mean expression per probe across all smaples
tmp_mat = cbind(exprs(annotated_probe_eset), mean_per_probe)

# per gene picking the probe with the highest mean
gene_mat = applyBy(tmp_mat,BY =  probe_factor, MARGIN = 2L, FUN =  function(mat, ...) 
{
  max_i = which.max(mat[,ncol(mat)])
  mat[max_i,]
})
gene_mat = gene_mat[,-ncol(gene_mat)]  # removing the mean_probe column

message("3) Computing Z-scores and FC")
dis_samp = grep("d_*", colnames(gene_mat))
cont_samp = grep("c_*", colnames(gene_mat))
FC_vec = apply(gene_mat, 1L, function(row) {mean(row[dis_samp], na.rm = T) - mean(row[cont_samp], na.rm = T)})
Z_gene_mat = apply(gene_mat, 2L, function(col) {(col - mean(col, na.rm = T))/sd(col, na.rm = T) })

message("4) Computing Z-test")
dis_n = length(dis_samp)
con_n = length(cont_samp)
cont_sd = apply(Z_gene_mat[,cont_samp],1L, sd, na.rm = T)
dis_sd = apply(Z_gene_mat[,dis_samp],1L, sd, na.rm = T)

n= nrow(Z_gene_mat)
Z_test = sapply(X = rownames(Z_gene_mat),USE.NAMES = F, FUN =  function(g)
{
  numerator = mean(Z_gene_mat[g, dis_samp], na.rm = T) - mean(Z_gene_mat[g, cont_samp], na.rm = T)
  denominator = sqrt(((cont_sd[g]^2)/con_n) + ((dis_sd[g]^2)/dis_n))
  numerator/denominator
})

message("5) Standardizing Z-test")
Ztest_mean = mean(Z_test)
Ztest_mean_sd = sd(Z_test)
Z_test_standard = (Z_test - Ztest_mean)/Ztest_mean_sd

message("6) Merging FC and Ztest and ploting")
comb_data = merge(FC_vec, Z_test_standard, by=0) 
colnames(comb_data)=c("gene", "FC", "Ztest")
rownames(comb_data) = comb_data[,"gene"]
comb_data = comb_data[,-1]
# hist(comb_data[,"FC"],100, main="FC")
# hist(comb_data[,"Ztest"],100, main="Ztest")
# plot(comb_data[,"FC"], comb_data[,"Ztest"], )
# comb_data[1:13,]
setwd(paste0("/storage/md_shenorr/inbaltz/microarrays/",RFST,"/",disease,"/"))
filename=paste0(gsub("data.rds", "", RDSfile),"all_orthologs.rds")
saveRDS(comb_data, filename) 
