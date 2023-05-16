# View Results
View a browsable version of the Protein-Metabolite associations as a Shiny app at: https://mbenson.shinyapps.io/protein-metabolite

---
---

# Protein-Metabolite Association Studies Identify Novel Proteomic Determinants of Metabolite Levels in Human Plasma

Mark D. Benson<sup>1*</sup>, Aaron S. Eisman<sup>1, 2*</sup>, Usman A. Tahir<sup>1</sup>, Daniel H. Katz<sup>1</sup>, Shuliang Deng<sup>1</sup>, Debby Ngo<sup>1</sup>, Jeremy M. Robbins<sup>1</sup>, Alissa Hofmann<sup>1</sup>, Xu Shi<sup>1</sup>, Shuning Zheng<sup>1</sup>, Michelle Keyes<sup>1</sup>, Zhi Yu<sup>3</sup>, Yan Gao<sup>4</sup>, Laurie Farrell<sup>1</sup>, Dongxiao Shen<sup>1</sup>, Zsu-Zsu Chen<sup>1</sup>, Daniel E. Cruz<sup>1</sup>, Mario Sims<sup>4</sup>, Adolfo Correa<sup>4</sup>, Russell P. Tracy<sup>5</sup>, Peter Durda<sup>5</sup>, Kent D. Taylor<sup>6</sup>, Yongmei Liu<sup>7</sup>, W. Craig Johnson<sup>8</sup>, Xiuqing Guo<sup>6</sup>, Jie Yao<sup>6</sup>, Yii-Der Ida Chen<sup>6</sup>, Ani W. Manichaikul<sup>9,10</sup>, Deepti Jain<sup>11</sup>, Qiong Yang<sup>12</sup>, NHLBI Trans-Omics for Precision Medicine (TOPMed) Consortium, Claude Bouchard<sup>13</sup>, Mark A. Sarzynski<sup>14</sup>, Stephen S. Rich<sup>9</sup>, Jerome I. Rotter<sup>6</sup>, Thomas J. Wang<sup>15</sup>, James G. Wilson<sup>1</sup>, Clary B. Clish<sup>3</sup>, Indra Neil Sarkar<sup>2</sup>, Pradeep Natarajan<sup>3</sup>, 16, 17</sup>, and Robert E. Gerszten<sup>1, 3†</sup>

1.	Division of Cardiovascular Medicine, Beth Israel Deaconess Medical Center, Boston, MA; 
2.	Center for Biomedical Informatics, Brown University, Providence, RI; 
3.	Broad Institute of Harvard and MIT, Cambridge, MA; 
4.	University of Mississippi Medical Center, Jackson, MS; 
5.	Department of Pathology Laboratory Medicine, Larner College of Medicine, University of Vermont, Burlington, VT; 
6.	The Institute for Translational Genomics and Population Sciences, Department of Pediatrics, The Lundquist Institute for Biomedical Innovation at Harbor-UCLA Medical Center, Torrance, CA; 
7.	Department of Medicine, Division of Cardiology, Duke Molecular Physiology Institute, Duke University Medical Center, Durham, NC; 
8.	Department of Biostatistics, University of Washington, Seattle, WA; 
9.	Center for Public Health Genomics, University of Virginia, Charlottesville, VA; 
10.	Division of Biostatistics and Epidemiology, Department of Public Health Sciences, University of Virginia, Charlottesville, VA; 
11.	University of Washington, Seattle, WA; 
12.	Department of Biostatistics, Boston University School of Public Health, Boston, MA;
13.	Human Genomic Laboratory, Pennington Biomedical Research Center, Baton Rouge, LA; 
14.	Department of Exercise Science, University of South Carolina, Columbia, SC; 
15.	Department of Medicine, UT Southwestern Medical Center, Dallas, TX; 
16.	Cardiovascular Research Center, Massachusetts General Hospital, Boston, MA; 
17.	Department of Medicine Harvard Medical School, Boston, MA

*indicates equal contribution in work and shared first co-authorship
†corresponding author

## Paper Abstract
While many novel gene-metabolite and gene-protein associations have been identified using high throughput biochemical profiling, systematic studies that leverage human genetics to illuminate novel protein-metabolite associations are lacking.  Here, we performed protein-metabolite association studies in 3,626 plasma samples from three population studies.  We detected over 171,800 protein-metabolite pairwise correlations between 1265 proteins and 365 metabolites (q-value ≤ 0.05).  Our findings include correlations from central metabolic and signaling pathways, such as the protein thyroxine binding globulin and the metabolite thyroxine (r = 0.51, q-value ≤ 1.0 x 10-300), as well as thousands of novel findings.  In enrichment analyses, we detected associations between plasma proteins and specific classes of amino acid, lipid and nucleotide metabolites, including associations between lipids and cathepsin proteases, serpin peptidase inhibitors, and secreted glycoproteins.  In Mendelian Randomization analyses, we identified putative causal protein-to-metabolite associations, including associations between the scavenger receptor CD36 and circulating levels of dozens of lipid species such as arachidonic acid (beta 0.14, q-value 2.0 x 10-4) and docosahexaenoic acid (beta 0.16, q-value 1.1 x 10-4).  We experimentally validated many of the top protein-to-metabolite associations in proof-of-concept plasma metabolomics studies in three murine knockout models.  These analyses identified previously unrecognized associations between bioactive proteins and metabolites in human plasma.  We provide publicly available data to be leveraged for studies in human metabolism and disease.

## About this repository
This repository contains the code used to generate results for the above citation organized into analysis steps by folders.

### data_prep
Protein, metabolite, and phenotype data from the Jackson Heart Study (JHS), the Multi-Ethnic Study of Atherosclerosis (MESA), and the Health, Risk Factors, Exercise Training and Genetics (HERITAGE) Family Study were ingested and put into a standardized form using the ``load_transform'' scripts within the data_prep folder. In addition, a key was generated to translate between metabolite names across the three studies. Finally, a key was generated to translate between protein IDs in HERITAGE and the other two studies (JHS and MESA).

### analysis
The analysis of these data was broken up into three parts. The first was a correlation analysis using the above-mentioned standardized and transformed versions of the raw protein and metabolite measurements. Next, a series of scripts were written to identify and extract genetic associations from previously performed protein and metabolite GWAS results that had been conducted on each of the three studies for mendelian randomization. Finally, those mendelian randomization analyses were meta-analyzed across studies.

### tables
Scripts to generate all tables, including supplemental tables, are provided.

### figures
Scripts to generate core components of figures are provided. Additional stylistic modifications (e.g. arrangement and font sizes) using Adobe Indesign and are not included here.

### shinyapp
A Shiny data browser and associated data files are available for download. Fully contained within the shinyapp folder is the R code and data necessary to browse the detailed results on your personal computer.
