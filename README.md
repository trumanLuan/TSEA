# 1. TSEA introduction
Genome-wide association studies (GWAS) and next-generation sequencing technologies have identified hundreds of thousands of disease-associated variants and genes. Interpretation of these variants could be greatly enhanced in tissue-specific systems. However, there are many diseases or traits where the causal tissues or cell types remain unknown. Tissue-Specific Enrichment Analysis (TSEA) is an R package to identify the most relevant tissues for candidate genes or for gene expression profiles. TSEA builds on two pre-processed reference panels. We implemented different statistic tests for different forms of query data. 
# 2. Usage
## 2.1 Installing TSEA
### Requirements
TSEA relies on R (>= 3.4), pheatmap (>= 1.0.10), RColorBrewer (>= 1.1)  
The pheatmap relies on CRAN. Please follow their installation instruction.  
`> install.packages("pheatmap")  `
### To download the codes, please do:
`git clone https://github.com/bsml320/TSEA.git`  
`cd TSEA`  
Then open the R:   
`R`  
`> install.packages("TSEA_1.0.tar.gz")  `
### TSEA loading
load the TSEA package and dependent library  
`> library(TSEA)`  
`> library(pheatmap)`  
## 2.2 Built-in data loading
load the t-statistic matrix for the GTEx panel  
`> load("data/GTEx_t_score.rda")`  
load the z-score matrix for the ENCODE panel  
`> load("data/ENCODE_z_score.rda")`  
Then "GTEx_t_score" and "ENCODE_z_score" panels will be loaded to R enviroment.  
## 2.3 Input data
TSEA deals with two types of enrichment analysis for different forms of query data. For convenience, we provide two TSEA functions for query gene lists (single sample and multiple sample), and another function for RNA-Seq expression profiles.    
### 2.3.1 TSEA for gene lists
When the query data are lists of genes, the Fisher’s Exact Test is implemented. The function is tsea.analysis(). The input is a vector of gene symbols. Here we used disease-associated genes identified from GWAS summary statistics as an example. The gene symbols can be found here:  
Load gene symbol from TSEA package.  
`> load("data/GWAS_gene.rda")`  
`> query.genes = GWAS_gene`  
Or you can read your own gene symbol list from a text file.  
`> dat = read.table("data/Gene_list.txt",head = F)`  
`> query.genes = dat[,1]`  
Nextly, we perform tissue-specific enrichment analysis for query gene list.  
`> tsea_t = tsea.analysis(query.genes, GTEx_t_score, ratio = 0.05, p.adjust.method = "bonferroni")`  
Here, the ratio is a value to define tissue-specific genes (default is 5%) and provides the first way of categorizing genes.  
The second way of grouping genes is based on the query genes. The two ways of category form a two by two table, which is used in the Fisher’s Exact Text.  
The Fisher's Exact Test results between query gene list and each tissue specific genes will be stored in variable `tsea_t`.  
You can check tissue-specific enrichment analysis result by:    
`> head(tsea_t)`  
For better visualization and summary, we provide one plot and one summary function to list the top 3 enriched tissues, simply run:  
`> tsea.plot(tsea_t, 0.05)`  
`> tsea_t_summary = tsea.summary(tsea_t)`  

### 2.3.2 TSEA for multiple gene lists  
In most condition, you might want to analysis multiple samples together, then you can upload a 0~1 table. In the table, gene labeled with 1 indicated significant associate within a sample, while 0 indicated not in a given sample. You can check the format of example data.  
Load multiple gene symbol from TSEA package.  
`> load("data/GWAS_gene_multiple.rda")`  
`> query.gene.list = GWAS_gene_multiple`  
Or you can read your own gene symbol list from a text file.  
`> dat = read.table("Gene_list_multiple.txt", head = T, row.names = 1)  `
`> query.gene.list = dat`  
Chech the total genes number for each sample. To keep result reliable, please keep at least 20 genes for each samples. 
`> colSums(query.gene.list)`  
Then, we can make tissue specific enrichment analysis for multiple samples by `tsea.analysis.multiple()` and plot the result by `tsea.plot()`. You can summary the top 3 most associated tissues by `tsea.summary() function` and save your result in to a text-format spreadsheet:  
Tissue-specific enrichment analysis in GTEx panel:  
`> tsea_t_multi = tsea.analysis.multiple(query.gene.list, GTEx_t_score, 0.05, p.adjust.method = "BH")`  
Save tissue-specific enrichment analysis result:  
`> write.csv(tsea_t_multi,"GWAS_multi_TSEA_in_GTEx_panel.csv")`  
Save the tissue-specific enrichment analysis plot:  
`> pdf ("GWAS_multi_TSEA_in_GTEx_panel.pdf",6,6,onefile = FALSE)`  
`> tsea.plot(tsea_t_multi, 0.05)`
`> dev.off()`   
Save your result in to a spreadsheet:  
`> tsea_t_multi_summary = tsea.summary(tsea_t_multi)`  
`> write.csv(tsea_t_multi_summary,"GWAS_multi_summary_GTEx_panel.csv")`

### 2.3.3 TSEA for RNA-seq profiles
For a quick start, user can use ENCODE example RNA-seq profiles:  
Load ENCODE query data:  
`> load("data/query_ENCODE.rda")`  
`> query.matrix = query_ENCODE`  
Load correction variable:  
`> load("data/correction_factor.rda")` 
As RNA-Seq samples are often heterogeneous, before in-depth analysis, it’s necessary to decode tissue heterogeneity to avoid samples with confounding effects. However, the raw discrete RPKM value should be normalized to continuous variable meet the normal distribution before t-test. We provided two normalization approaches: `"z-score"` and `"abundance"` in function `tsea.expression.normalization()`:  
(1) `z-score` normalization will calculate a z-score for the query sample for each tissue in the reference panel as below: e_i=(e_0-μ_t))/sd_t, where μ_t and sd_t were the mean and SD of tissue t.   
(2) `abundance` normalization will provide an abundance correction approach for the query sample for each tissue in the reference panel as below: e_i=(log2(e_0+1)/(log2(u_t+1)+1).  
We have the preloaded the test RPKM variable in `query.matrix` and correction variable in `correction_factor`, we take "abundance" normalization approach as an example, simply type:  
RNA-Seq profiles scale by abundance normalization:
`> query_mat_abundance_nor = tsea.expression.normalization(query.matrix, correction_factor, normalization = "abundance")`  
After get normalized RPKM value, we submit it for `tsea.expression.decode()`:  
`> tseaed_in_GTEx = tsea.expression.decode(query_mat_abundance_nor, GTEx_t_score, 0.05, p.adjust.method = "BH")`  
Then, the tissue specific enrichment analysis for query RNA-seq is finish. After tissue specific enrichment decode analysis, one-side t-test results between query RNA-seq sample tissue specific genes (top 5%) versus remains genes (95%) is stored in variable `tseaed_in_GTEx`. Further analysis for top 3 most associated tissues is similar to previous analysis:  
`> tsea.plot(tseaed_in_GTEx, 0.05)`  
`> tseaed_in_GTEx_summary = tsea.summary(tseaed_in_GTEx)`  
`> write.csv(tseaed_in_GTEx_summary,"RNAseq_summary_in_GTEx_panel.csv")`  

## Citation
Pei G., Dai Y., Zhao Z, Jia P. (2018) Tissue-Specific Enrichment Analysis (TSEA) to decode tissue heterogeneity. In submission.  


















