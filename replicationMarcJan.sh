#!/bin/bash

ml Java


#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name="InteractionAnalysis"

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
#ml nextflow

# We set the following variables for nextflow to prevent writing to your home directory (and potentially filling it completely)
# Feel free to change these as you wish.
export SINGULARITY_CACHEDIR=../singularitycache
export NXF_HOME=../nextflowcache

# Disable pathname expansion. Nextflow handles pathname expansion by itself.
set -f

c="LL"


# Define paths
base_folder=/groups/umcg-bios/tmp04/projects/BIOS_for_eQTLGenII/pipeline/20220426/
# Genotype data
# Full path to the folder with imputed filtered vcf files produced by eQTLGen pipeline 2_Imputation step (postimpute folder)]
vcf_dir_path=${base_folder}/2_Imputation/out/${c}/postimpute/

# raw expression data (same as input to DataQC step)
raw_exp_path=/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/data/${c}/${c}_raw_expression.txt.gz
# normalized expression data (output of the DataQC step)
norm_exp_path=${base_folder}/1_DataQC/out/${c}/outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt
# File that contains cohort covariates: E.g. sex and age.
covariate_path=/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/BIOS_covariates.txt
# genotype to expression coupling file
gte_path=/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/output/2023-03-16-sex-specific-analyses/run1/data/${c}/${c}.gte

# covariate to test (name of the sex column in the covariate file)  gender_F1M2
covariate_to_test=gender

# Path to genotype PCs (output of dataQC step)
genotype_pcs_path=${base_folder}/1_DataQC/out/${c}/outputfolder_gen/gen_PCs/GenotypePCs.txt
# Path to expression PCs (output of dataQC step)
expression_pcs_path=${base_folder}/1_DataQC/out/${c}/outputfolder_exp/exp_PCs/exp_PCs.txt

# output folder
output_path=/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/replicationMarcJan1/


# Path to the nextflow interaction analysis folder
script_folder=/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/eQTLGenInteractions/

qtls_to_test=/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/sig_Qtls_sceQTLGen_062024_3.txt.gz
chunk_file=${script_folder}/data/ChunkingFile.txt
exp_platform=RNAseq #options: RNAseq; RNAseq_HGNC; HT12v3; HT12v4; HuRef8; AffyU219; AffyHumanExon

expression_eigenvectors=/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/downloadData/EigenvectorsTop1000.txt.gz # The expression eigenvectors as calculated using all eqtlgen samples
expression_ics=/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/downloadData/Ica100.txt.gz # The expression independent components as calculated using all eqtlgen samples


mkdir -p ${output_path}

# Command:
NXF_VER=24.04.4 ../nextflow/nextflow-24.04.4-all run /groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/eQTLGenInteractions/InteractionAnalysis.nf \
--vcf_dir $vcf_dir_path \
--raw_expfile ${raw_exp_path} \
--norm_expfile ${norm_exp_path} \
--gte ${gte_path} \
--covariates $covariate_path \
--exp_platform ${exp_platform} \
--cohort_name ${cohort_name} \
--covariate_to_test $covariate_to_test \
--qtls_to_test $qtls_to_test \
--genotype_pcs $genotype_pcs_path \
--expression_eigenvectors $expression_eigenvectors \
--expression_ics $expression_ics \
--chunk_file $chunk_file \
--outdir ${output_path}  \
--run_stratified false \
--preadjust false \
--cell_perc_interactions false \
-profile slurm \
--dev false \
-resume




# -resume \
#-profile slurm
#standard
