# Soy-Meta-GWAS
Soybean Meta-GWAS related files
TABLE OF CONTENTS
1) INTRODUCTION
2) MAIN FILES
3) SUPPLEMENTAL FILES

1. INTRODUCTION
This page contains related to the analysis of NPGS germplasm characterization studies, as well as the meta-analysis of these studies.

2. FILES
Phenotypic data in numeric format is compiled in "Phenotype data.csv". The first column, "acid" is the accession id, colloquially known as the PI number, for each soybean accession included in one or more studies. Each subsequent column is identified by the abbreviated study name generated by NPGS, as well as the trait to which the data pertains. 

Some of the initial Technical Bulletins that were used to generate the phenotypic data are provided, and can be identified by the naming scheme "Tech####.pdf", where the #### is replaced by the Technical Bulletin number. These files are provided here primarily due to the authors' difficulty in locating these files initially.

"VCR Test Definition.xlsx" provides a list of all studies included in the corresponding paper, as well as the experiment number from NPGS (ENO), name from NPGS (ENAME), methods used (METHODS), literature id (LITID), and a binary flag to show inclusion in this Github page.

In addition, four R scripts are provided related to the original GWAS analysis (non-meta-analysis portion). "Preparation_GD_GM_from_soybase_vcf.R" was used to create the master file for genotypic data used in each study. "GWAS prep.R" was used to convert the imputed master file to a format compatible with GAPIT. Several tests were performed with different model types to determine which model appeared to give the most accurate set of results, before moving that model to "NPGS GWAS.R" for complete analysis. Finally, "CompileSignificantSNPs.R" was used to identify associations which were significant at the study-specific Bonferroni threshold.

3. SUPPLEMENTAL FILES
Supplemental files from the submission are provided for convenience here as well. 

Supplemental Figure 1 shows LD within the 302 resequencing panel (Zhou et al. 2016) for Dt1, PMDH1, and FATB1a. LD is calculated from the lead SNP identified by the largest study which detected each locus.

Supplemental Figure 2 shows LD decay surrounding the three previously listed genes based on multiple studies which detected the association, using SoySNP50k chip marker data for that study.

Supplemental Figure 3a shows all significant SNPs detected in individual studies for seed composition traits (outer ring), agronomic traits (middle ring), and disease traits (inner ring).

Supplemental Figure 3b shows all significant SNPs detected in meta-GWAS for seed composition traits (outer ring), agronomic traits (middle ring), and disease traits (inner ring).

Supplemental Table 1 contains all significant SNPs (Based on Bonferroni criteria for individual studies, or implied Bonferroni level of 0.05/42080 for meta-GWAS) found within this study, as well as the P-value, minor allele frequency, number of observations (nobs), false discovery rate adjusted P-value, candidate gene, QTL number (this study), and other studies that have detected this association (limited).

Supplemental Table 2 contains the evaluation number (ENO), evaluation name (ENAME), provided methods (NPGS), reference, and a list of traits studies within each experiment.

Supplemental Table 3 defines traits and lists the number of unique loci and candidate genes identified for each trait.

Supplemental Table 4 illustrates the differences in map position for markers associated with the SACPD-C gene conferring altered stearic acid content between Wms82.a1 and Wms82.a2.

Supplemental Table 5 summarizes the lead SNP for each experiment/trait/locus combination, the candidate gene identified for the locus, whether the locus was found in only the Meta-GWAS, only in individual studies, or can be found in both, and the D' value for the lead SNP from that study with the candidate gene based on SoySNP50k chip marker data (if only Comparison SNP column given, a marker within the candidate gene, if Comparison SNP and Comparison SNP2 both given, immediate flanking markers for that gene).