# Information about the repository
# 

This directory contains all computer code to perform the analyses
from our upcoming manuscript, "A polygenic score-based approach
to identify gene-drug interactions stratifying breast cancer risk".

**June 14, 2021: Please note that the scripts could use further commenting and documentation. This is on the to-do list. Currently, this mostly serves as a repository for myself to access any relevant scripts but other researchers can access these scripts in their current format as well.**

All scripts are written in the R or Bash programming languages. The analysis was performed on a linux system. Plots were created on a macOS Catalina system.

If you can not find the code you are looking for or have any questions, please contact:

Andrew Marderstein
anm2868@med.cornell.edu

If you use our work, please cite:

Marderstein, A.R., Kulm, S., Peng, C., Tamimi, R.M., Clark, A.G.^, Elemento, O.^ (2021). A polygenic score-based approach to identify gene-drug interactions stratifying breast cancer risk. *medRxiv*. doi: https://doi.org/10.1101/2021.05.03.21256511.

#### Generating the analysis data

The create_data directory contains data-generating files.

#### Calculating polygenic scores

Please see the create_score and create_score_PGS_catalog directories.

#### Polygenic score metrics

	pgs_select_optimal_score_plot_plus_heatmap.R
	pgs_corr_plot.R
	identify_best_score.R
	cv_overfit.R
	interaction_analysis_plots.R

#### Polygenic-drug interaction tests

	pgs_drug_gxdrug.R
	pgs_qq.R
	variation_explained_by_gxdrug.R
	pgs_variation_explained.R


#### Single-SNP interaction analyses

Creating genotype files for the single-SNP analysis using:

	create_genotype_files_for_snp_gxe.sh
	snp_drug_gxdrug.R
	snp_drug_gxdrug.analyze_res.R
	coord_int_statistics.R
	create_snpxcort_supp_tab.R
	survival_plots.R

#### In silico single-SNP functional analyses:

	ukb-id_to_open-targets-id.sh
	enrichment_analysis.R
	enrichr_analyses.R
	gtex_snp.R
	
#### Sensitivity analyses for confounders:

	confounding_analyses.R
	corticosteroid_indic_create.R
	sensitivity_remove-indic.R

#### Other scripts

	figure1.R
