# 1. TSEA introduction
Genome-wide association studies (GWAS) and next-generation sequencing technologies have identified hundreds of thousands of disease-associated variants and genes. Interpretation of these variants could be greatly enhanced in tissue-specific systems. However, there are many diseases or traits where the causal tissues or cell types remain unknown. Tissue-Specific Enrichment Analysis (TSEA) is an R package to identify the most relevant tissues for candidate genes or for gene expression profiles. TSEA builds on two pre-processed reference panels. We implemented different statistic tests for different forms of query data. 
# 2. Usage
## 2.1 Installing TSEA
### Requirements
TSEA relies on R (>= 3.4), pheatmap (>= 1.0.10), RColorBrewer (>= 1.1)
The pheatmap relies on CRAN. Please follow their installation instruction.
install.packages("pheatmap")
### To download the codes, please do:
git clone https://github.com/bsml320/TSEA.git
cd TSEA
Then open the R：
R
> install.packages("TSEA_1.0.tar.gz")
### TSEA loading
load the TSEA package and dependent library
> library(TSEA)
> library(pheatmap)
## 2.2 Built-in data loading
load the t-statistic matrix for the GTEx panel, 
> load("data/GTEx_t_score.rda")
then "GTEx_t_score" will be loaded to R enviroment

load the z-score matrix for the ENCODE panel
> load("data/ENCODE_z_score.rda")
then "ENCODE_z_score" will be loaded to R enviroment
## 2.3 Input data
TSEA deals with two types of enrichment analysis for different forms of query data. 
### 2.3.1 TSEA for gene lists
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

### 2.3.2 TSEA for multiple gene lists
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

### 2.3.3 TSEA for RNA-seq profiles
For a quick start, user can use ENCODE example RNA-seq profiles:
Load ENCODE query data
> data(query_ENCODE)
> query.matrix = query_ENCODE
> head(query.matrix)[,1:4]

Adrenal Gland Body of Pancreas Breast Epithelium Camera-type Eye 
TSPAN6	11.639167	5.3900	11.038333	24.6475
TNMD		0.010000	0.1475	2.243333	12.4325
DPM1		18.819167	9.8125	14.215000	24.0250
SCYL3		3.812500	2.5925	5.626667	10.4975
C1orf112	1.643333	2.3075	2.858333	14.6125
FGR		5.282500	1.6750	8.476667	1.5625

As RNA-Seq samples are often heterogeneous, before in-depth analysis, it’s necessary to decode tissue heterogeneity to avoid samples with confounding effects. However, the raw discrete RPKM value should be normalized to continuous variable meet the normal distribution before t-test. We provided two normalization approaches: "z-score" and "abundance" in function tsea.expression.normalization():
	z-score normalization will calculate a z-score for the query sample for each tissue in the reference panel as below: e_i=(e_0-μ_t))/sd_t, where μ_t and sd_t were the mean and SD of tissue t. 
	abundance normalization will provide an abundance correction approach for the query sample for each tissue in the reference panel as below: e_i=(log2(e_0+1)/(log2(u_t+1)+1).

We have the preloaded the test RPKM variable in query.matrix and correction variable in correction_factor, we take "abundance" normalization approach as an example, simply type:
RNA-Seq profiles scale by abundance normalization
query_mat_abundance_nor = tsea.expression.normalization(query.matrix, correction_factor, normalization = "abundance")
> head(query_mat_abundance_nor)[,1:4]
Adrenal Gland Body of Pancreas Breast Epithelium Camera-type Eye
C1orf112	0.7427561	0.9140352	1.0317423 	2.0998577
FGR 		0.4884476	0.2615171	0.5977016 	0.2500974
CFH    	0.8805367	0.5989594 	0.8540355 	0.8547878
FUCA2  	1.0366443	0.9182102	0.7530246	0.7801565
NFYA		0.6386405	0.5508964	0.8167324	1.1198416
STPG1		0.9974109	0.3903567	0.5879190	0.8472181

After get normalized RPKM value, we submit it for tsea.expression.decode()
> tseaed_in_GTEx = tsea.expression.decode(query_mat_abundance_nor, 
		GTEx_t_score, 0.05, p.adjust.method = "BH")
> head(tseaed_in_GTEx)[,1:3]
            Adrenal Gland Body of Pancreas Breast Epithelium
Adipose - Subcutaneous	7.093272e-49 4.636686e-44 5.779239e-142
Adipose - Visceral (Omentum) 2.199051e-33	9.641733e-32 1.313532e-112
Adrenal Gland	9.925492e-220 3.404965e-33 1.476023e-25
Artery - Aorta	2.749910e-39 1.075081e-21 5.919623e-47
Artery - Coronary	6.702724e-42 4.048237e-22 6.188863e-44
Artery - Tibial	1.461867e-37 7.482580e-24 5.141699e-52 

Then, the tissue specific enrichment analysis for query RNA-seq is finish. After tissue specific enrichment decode analysis, one-side t-test results between query RNA-seq sample tissue specific genes (top 5%) versus remains genes (95%) is stored in variable tseaed_in_GTEx. Further analysis for top 3 most associated tissues is similar to previous analysis

> tsea.plot(tseaed_in_GTEx, 0.05)
> tseaed_in_GTEx_summary = tsea.summary(tseaed_in_GTEx)
> write.csv(tseaed_in_GTEx_summary,"RNAseq_summary_in_GTEx_panel.csv")

## Citation
Pei G., Dai Y., Zhao Z, Jia P. (2018) Tissue-Specific Enrichment Analysis (TSEA) to decode tissue heterogeneity. Bioinformatics, in submission.


















