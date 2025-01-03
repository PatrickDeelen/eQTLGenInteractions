#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage() {
  log.info"""
  =======================================================
    Interaction analysis v${workflow.manifest.version}
  =======================================================
    
  """.stripIndent()
}

/*
 * Parameters
 */

params.bfile = ''
params.vcf_dir = ''
params.bgen_dir = ''

params.genes_to_test = ''
params.qtls_to_test = ''
params.lab_cell_perc = ''

params.signature_matrix_name = "LM22"
params.deconvolution_method = "dtangle"
params.num_perm = 1

params.run_stratified = false
params.preadjust = false
params.cell_perc_interactions = false
params.expr_pcs = ''
params.num_expr_PCs = 20

//params.expressionEigenvectors = '' // The expression eigenvectors as calculated using all eqtlgen samples
//  params.expressionIcs = '' // The expression independent components as calculated using all eqtlgen samples

/*
 * Channel declarations
 */
raw_expr_ch = Channel.fromPath(params.raw_expfile)
filt_exp_ch = Channel.fromPath(params.norm_expfile) // the expression data normalized as used for the eqtlgen eqtl mapping
outdir_ch = Channel.fromPath(params.outdir, type: 'dir')
covars_ch = Channel.fromPath(params.covariates)
gte_ch = Channel.fromPath(params.gte)

annotation_ch = Channel.fromPath("$projectDir/data/LimixAnnotationFile.GRCh38.110.txt.gz")
gene_lengths_ch = Channel.fromPath("$projectDir/data/GeneLengths_GRCh38.110_ensg.txt.gz")

// chunk channel for the interaction analysis: chrom -> chunk
Channel
    .fromPath(params.chunk_file)
    .splitCsv( header: false )
    .map { row -> tuple(row[0].split(':')[0], row[0]) }
    .set { chunk_ch }

// genotype channel declaration. Interaction analysis doesn't run with bgen files yet
if (params.bgen_dir != '') {
  Channel.from(1..22)
  .map { chr -> tuple("$chr", file("${params.bgen_dir}/*chr${chr}.bgen"), file("${params.bgen_dir}/*chr${chr}.sample")) }
  .ifEmpty { exit 1, "Input .bgen files not found!" }
  .set { chr_bgen_pairs }
} else if (params.vcf_dir != '') {
  Channel.from(1..22)
  .map { chr -> tuple("$chr", file("${params.vcf_dir}/*chr${chr}.filtered.vcf.gz")) }
  .ifEmpty { exit 1, "Input .vcf.gz files not found!" }
  .set { chr_vcf_pairs }
  Channel.empty()
    .set { chr_bgen_pairs }
} else {
  Channel
    .from(params.bfile)
    .ifEmpty { exit 1, "Input plink prefix not found!" }
    .map { genotypes -> [file("${genotypes}.bed"), file("${genotypes}.bim"), file("${genotypes}.fam")]}
    .set { bfile_ch }
} 


include { PREPARE_COVARIATES; NormalizeExpression; ConvertVcfToBgen; MergeBgenPerChr; MapEigenvectorsAndIcs; Transpose } from './modules/prepare_data.nf'
include { RUN_INTERACTION_QTL_MAPPING; IeQTLmapping } from './modules/interaction_analysis3.nf'
include { RUN_STRATIFIED_ANALYSIS; RunEqtlMappingPerGenePlink } from './modules/stratified_analysis.nf'

/* 
 * Analysis
 */

workflow {
    params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }


    if(params.dev){
        // this block is only run if test flag is used for development
        NormalizeExpression(raw_expr_ch, filt_exp_ch, params.exp_platform, gte_ch )
        norm_exp_ch = NormalizeExpression.out.norm_expression_table
        covariates_ch = PREPARE_COVARIATES(params.exp_platform, raw_expr_ch, norm_exp_ch, params.signature_matrix_name, params.deconvolution_method,covars_ch, gene_lengths_ch, annotation_ch, Channel.fromPath(params.genotype_pcs), Channel.fromPath(params.gte), filt_exp_ch)


    } else {


        /*
         * Prepare normalized expression and covariate data
         */
        NormalizeExpression(raw_expr_ch, filt_exp_ch, params.exp_platform, gte_ch )
        norm_exp_ch = NormalizeExpression.out.norm_expression_table
        covariates_ch = PREPARE_COVARIATES(params.exp_platform, raw_expr_ch, norm_exp_ch, params.signature_matrix_name, params.deconvolution_method,covars_ch, gene_lengths_ch, annotation_ch, Channel.fromPath(params.genotype_pcs), Channel.fromPath(params.gte), filt_exp_ch)


        RUN_INTERACTION_QTL_MAPPING(Transpose(filt_exp_ch   , "expressionT.txt"), covariates_ch, annotation_ch, Channel.of(params.covariate_to_test), chunk_ch.map { it[1] }, Channel.fromPath(params.qtls_to_test),  ConvertVcfToBgen(params.vcf_dir).bgen_ch)



    }
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
