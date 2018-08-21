# 1. TSEA introduction
Genome-wide association studies (GWAS) and next-generation sequencing technologies have identified hundreds of thousands of disease-associated variants and genes. Interpretation of these variants could be greatly enhanced in tissue-specific systems. However, there are many diseases or traits where the causal tissues or cell types remain unknown. Tissue-Specific Enrichment Analysis (TSEA) is an R package to identify the most relevant tissues for candidate genes or for gene expression profiles. TSEA builds on two pre-processed reference panels. We implemented different statistic tests for different forms of query data. 
# 2. Usage
# 2.1 Installing TSEA
## Requirements
TSEA relies on R (>= 3.4), pheatmap (>= 1.0.10), RColorBrewer (>= 1.1)
The pheatmap relies on CRAN. Please follow their installation instruction.
install.packages("pheatmap")
## To download the codes, please do:
git clone https://github.com/bsml320/TSEA/TSEA_1.0.tar.gz
install.packages("TSEA")
## TSEA loading
load the TSEA package and dependent library
> library(TSEA)
> library(pheatmap)
# 2.2 Built-in data loading
load the t-statistic matrix for the GTEx panel
> data(GTEx_t_score)
load the z-score matrix for the ENCODE panel
> data(ENCODE_z_score)
# 2.3 Input data
TSEA deals with two types of enrichment analysis for different forms of query data. 
# 2.3.1 TSEA for gene lists
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
> tsea.plot(tsea_t, 0.05)
> tsea_t_summary = tsea.summary(tsea_t)
> head (tsea_t_summary)
Top1	tissue1_p-value	Top2	tissue2_p-value	Top3	tissue3_p-value                              query	"Muscle - Skeletal"	"0.000213565633284591"	"Artery - Tibial" "0.00257812976565265"	"Adipose - Visceral (Omentum)"	"0.0109584998307628"

# 2.3.2 TSEA for multiple gene lists
In most condition, you might want to analysis multiple samples together, then you can upload a 0~1 table. In the table, gene labeled with 1 indicated significant associate within a sample, while 0 indicated not in a given sample. You can check the format of example data.
Load multiple gene symbol from TSEA package.
> data(GWAS_gene_multiple)
> query.gene.list = GWAS_gene_multiple
 Or you can read multiple gene symbol from a text file.
> dat = read.table("Gene_list_multiple.txt", head = T, row.names = 1)
> query.gene.list = dat
> query.gene.list [c(10:15),c(1:10)]
       ALZ ADHD ASD BD MDD SCZ BMI FN-BMD LS-BMD EDU
A4GNT    0    0   0  0   0   0   0      0      0   0
AA06     0    0   0  0   0   0   0      0      0   0
AAAS     0    0   0  0   0   0   0      0      1   0
AACS     0    0   0  0   0   0   0      0      0   0
AACSP1   0    0   0  0   0   0   0      0      0   0
AADAC    0    0   0  0   0   0   0      0      0   0

> colSums(query.gene.list)
  ALZ   ADHD    ASD     BD    MDD    SCZ    BMI FN-BMD LS-BMD    EDU HEIGHT    WHR     CD    IBD     RA     UC 
   121     20     20     25     20     88    233    109    107    487   3311     98    336    480    213    275 
   AAM    CAD     FG     FI    HDL    LDL     TC     TG    T1D    T2D 
   270     71     69     26    368    305    329    223    329     43

Then, we can make tissue specific enrichment analysis for multiple samples by tsea.analysis.multiple() and plot the result by tsea.plot() as showed in Fig. 2. You can summary the top 3 most associated tissues by tsea.summary() function and save your result in to a text-format spreadsheet, simply type:
Tissue-specific enrichment analysis in GTEx panel
> tsea_t_multi = tsea.analysis.multiple(query.gene.list, 
		GTEx_t_score, 0.05, p.adjust.method = "BH")
> #write.csv(tsea_t_multi,"GWAS_multi_TSEA_in_GTEx_panel.csv")
Heatmap Plot for TSEA result
> #pdf ("GWAS_multi_TSEA_in_GTEx_panel.pdf",6,6,onefile = FALSE)
> tsea.plot(tsea_t_multi, 0.05)
> #dev.off()
Save your result in to a spreadsheet
> tsea_t_multi_summary = tsea.summary(tsea_t_multi)
> write.csv(tsea_t_multi_summary,"GWAS_multi_summary_GTEx_panel.csv")








