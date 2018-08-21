## TSEA
Genome-wide association studies (GWAS) and next-generation sequencing technologies have identified hundreds of thousands of disease-associated variants and genes. Interpretation of these variants could be greatly enhanced in tissue-specific systems. However, there are many diseases or traits where the causal tissues or cell types remain unknown. Tissue-Specific Enrichment Analysis (TSEA) is an R package to identify the most relevant tissues for candidate genes or for gene expression profiles. TSEA builds on two pre-processed reference panels. We implemented different statistic tests for different forms of query data. 
## Installing TSEA
# Requirements
TSEA relies on R (>= 3.4), pheatmap (>= 1.0.10), RColorBrewer (>= 1.1)
The pheatmap relies on CRAN. Please follow their installation instruction.
install.packages("pheatmap")
# To download the codes, please do:
git clone https://github.com/bsml320/TSEA/TSEA_1.0.tar.gz
install.packages("TSEA")
## Usage
# TSEA loading
load the TSEA package and dependent library
> library(TSEA)
> library(pheatmap)
# Built-in data loading
load the t-statistic matrix for the GTEx panel
> data(GTEx_t_score)
load the z-score matrix for the ENCODE panel
> data(ENCODE_z_score)
## Input data
# TSEA for gene lists
When the query data are lists of genes, the Fisher’s Exact Test is implemented. The function is tsea.analysis(). The input is a vector of gene symbols. Here we used disease-associated genes identified from GWAS summary statistics as an example. The gene symbols can be found here:
Load gene symbol from TSEA package.
> data(GWAS_gene)
> query.genes = GWAS_gene

Or you can read gene symbol from a text file.
> dat = read.table("Gene_list.txt",head = F)
> query.genes = dat[,1]
Tissue-specific analysis for query gene list.
tsea_t = tsea.analysis(query.genes, GTEx_t_score, ratio = 0.05, p.adjust.method = "bonferroni")
Here, the ratio is a value to define tissue-specific genes and provides the first way of categorizing genes. The second way of grouping genes is based on the query genes. The two ways of category form a two by two table, which is used in the Fisher’s Exact Text.
The Fisher's Exact Test results between query gene list and each tissue specific genes will be stored in variable tsea_t.

Check tissue-specific enrichment analysis result.
> head(tsea_t)
                                  query
Adipose - Subcutaneous       1.00000000
Adipose - Visceral (Omentum) 0.01095850
Adrenal Gland                1.00000000
Artery - Aorta               0.21208614
Artery - Coronary            0.01095850
Artery - Tibial              0.00257813

For better visualization and summary, we provide one plot and one summary function to list the top 3 enriched tissues, simply run:
# TSEA result plot and summary
> tsea.plot(tsea_t, 0.05)
> tsea_t_summary = tsea.summary(tsea_t)
> head (tsea_t_summary)
Top1	tissue1_p-value	Top2	tissue2_p-value	Top3	tissue3_p-value                              query	"Muscle - Skeletal"	"0.000213565633284591"	"Artery - Tibial" "0.00257812976565265"	"Adipose - Visceral (Omentum)"	"0.0109584998307628"



