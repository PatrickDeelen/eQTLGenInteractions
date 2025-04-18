#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * run TPM normalization on the raw expression file and rename gene ids to gene names.
 * This expression table will be used in cell type deconvolution
 */
process TPM {
    tag "TPM normalization"
    label "short"

    input:
      path raw_expression
      path gene_lengths


    output:
      path "*.TPM.txt.gz"

    script:
    """
    Rscript $projectDir/bin/normalize_TPM.R \
        ${raw_expression} \
        ${gene_lengths} \
        expression.TPM.txt 
        
    gzip -f expression.TPM.txt
    """
}

/*
 * run deconvolution using $deconvolution_method method with $signature_matrix signature matrix
 */
process Deconvolution {
    tag "deconvolution"
    label "medium1"

    publishDir params.outdir, mode: 'copy'
    
    input:
      path tpm_expression
      path signature_matrix
      val deconvolution_method
      val exptype

    output:
      path "cell_counts.txt", emit: cellCounts
      path "cell_counts.txt.md5"

    script:
    """
    Rscript $projectDir/bin/run_deconvolution.R \
        ${tpm_expression} \
        ${signature_matrix} \
	    ${deconvolution_method} \
        cell_counts.txt \
        ${exptype}

    md5sum cell_counts.txt > cell_counts.txt.md5


    """
}

/*
 * Combine major covariates with cell counts, genotype PCs and RNA-quality. Plot covariate distributions.
 */
process CombineCovariatesRNAqual {
    label "medium1"
    debug false

    publishDir params.outdir, mode: 'copy'

    input:
        path general_covariates
        path cell_counts
        path genotype_PCs
        path gte
        path rna_qual
        path eigenAndIc
        path filt_exp_ch

    output:
        path ("covariates.combined.txt"), emit: covariates_ch
        path ("covariates.combined.txt.md5")
        path ("availableCovariates.txt")
        path ("availableCovariates.txt.md5")

    script:

      """
      echo "test6"
      Rscript $projectDir/bin/combine_all_covariates.R -s ${general_covariates} -c ${cell_counts} -g ${genotype_PCs} -i ${gte} -o covariates.combined.txt -r ${rna_qual} -v ${eigenAndIc} -e ${filt_exp_ch}
       md5sum covariates.combined.txt > covariates.combined.txt.md5
       md5sum availableCovariates.txt > availableCovariates.txt.md5
      """


}

/*
 * Calculate RNA quality score for RNA-seq datasets as the Pearson correlation between each sample and the average expression over all samples
 */
process CalculateRNAQualityScore {
    label "short"
    
    input:
    path tmm_expression_data
    
    output:
	path ("RNA_quality.txt_CorrelationsWithAverageExpression.txt.gz"), emit: rnaquality_ch

    script:
    """
	python3 $projectDir/bin/correlate_samples_with_avg_gene_expression.py -ex ${tmm_expression_data} -op RNA_quality.txt -log2
    """

}

/*
 * Normalize raw expression: use TMM for RNA-seq data and quantile normalization for array data. 
 * Do not apply sample and gene filtering again, instead use the normalized output from eQTLGen DataQC pipeline, where these filters have been applied. 
 * We can't use that DataQC output directly because we don't need log2 and INT that was applied there.
 */
process NormalizeExpression {
    label "medium2"
    publishDir params.outdir, mode: 'copy'
    
    input:
        path(raw_expr)
        path(norm_expr)
        val(exp_platform)
        path(gte)

    output:
	    path ('outputfolder_exp'), emit: expression_folder 
	    path ('outputfolder_exp/exp_data_preprocessed.txt'), emit: norm_expression_table
    
    shell:
    '''
       if [[ !{exp_platform} == "HT12v3" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_IlluminaHT12v3.txt
        elif [[ !{exp_platform} == "HuRef8" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_IlluminaHuRef8.txt
        elif [[ !{exp_platform} == "HT12v4" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_IlluminaHT12v4.txt
        elif [[ !{exp_platform} == "RNAseq" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_RNAseq.txt
        elif [[ !{exp_platform} == "AffyU219" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_AffyU219.txt
        elif [[ !{exp_platform} == "AffyHumanExon" ]]; then
            probe_mapping_file=!{baseDir}/data/EmpiricalProbeMatching_AffyHumanExon.txt
        elif [[ !{exp_platform} == "RNAseq_HGNC" ]]; then
            probe_mapping_file=!{baseDir}/data/HgncToEnsemblProbeMatching.txt
        fi
                
        outdir=${PWD}/outputfolder_exp/

        Rscript !{baseDir}/bin/ProcessExpression_v2.R \
           -e !{raw_expr} \
           -n !{norm_expr} \
           -l !{gte} \
           -p !{exp_platform} \
           -m $probe_mapping_file \
           -o ${outdir}


        cd ${outdir}

        md5sum * > outputfolder_exp.md5

    '''
}

/*
* Create samples scores based on the eigenvectors and independent components
*/
process MapEigenvectorsAndIcs {
    label "medium2"
    debug false


    input:
        path normalized_expression_data
        path expression_eigenvectors
        path expression_ics

    output:
        path "mappedEigenvectorsAndIcs.txt"


    script:
        """
        echo "Map eigen and ic"
        Rscript $projectDir/bin/mapEigenvectorsAndIcs.R -e ${normalized_expression_data} -v ${expression_eigenvectors} -i ${expression_ics} -o mappedEigenvectorsAndIcs.txt
        """
}

/*
 * Split the combined covariate table into 2 according to the covariate of interest (E.g. sex or age). 
 * In case of a binary covariate the 2 values will be used for the split (E.g. males and females separately), 
 * in case of a quantitative covariate, samples falling in the first and last quartile will be written to the output files.
 */
process SplitCovariates {
    label "medium1"
    publishDir params.outdir, mode: 'copy'
    
    input:
        path normalized_expression_data
        path combined_covariates
        val covariate_name
    

    output:
        path "*txt"

    script:
    """
        Rscript $projectDir/bin/split_covariate_into_bins.R $combined_covariates $covariate_name $normalized_expression_data ./ 
    """

}

process Transpose {
    label "medium1"

    input:
        path inputMatrix
        val outputFileName

    output:
        path outputFileName

    script:
    """
        echo $inputMatrix
        echo $outputFileName
        Rscript $projectDir/bin/transpose.R $inputMatrix $outputFileName
    """

}

/*
 * Convert VCF files to bgen (not currently used)
 */
process ConvertVcfToBgen {
    label "medium2"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true, failOnError: true

    input:
        path vcfDir

    output:
        tuple path("merged.bgen"), path("merged.bgen.metafile"), emit: bgen_ch
        path("merged.sample")
    
    script:
    if (params.qtls_to_test == ''){
               """

                   HOME=${PWD}/home
                   mkdir -p ${HOME}

echo 2
          java -jar /groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/GenotypeHarmonizer-1.4.28-SNAPSHOT/GenotypeHarmonizer.jar -i ${vcfDir} -I VCF_FOLDER -O BGEN -o merged --mafFilter 0.01 --hweFilter 1e-06 --callRateFilter 0.95 --machR2Filter 0.4 --genotypeField GP
          python -c "from bgen_reader import read_bgen; bgen = read_bgen('merged.bgen', verbose=False)"
          """

      }
      else {
       """
          zcat  ${params.qtls_to_test} | cut -f 2  > snpsToTest.txt

    export HOME="./home"
    echo "test"
    echo \${PWD}
        echo \${HOME}


    mkdir -p \${HOME}


            java -jar /groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/GenotypeHarmonizer-1.4.28-SNAPSHOT/GenotypeHarmonizer.jar -i ${vcfDir} -I VCF_FOLDER -O BGEN -o merged --mafFilter 0.01 --hweFilter 1e-06 --callRateFilter 0.95 --machR2Filter 0.4 --genotypeField GP --variantFilterList snpsToTest.txt
            python -c "from bgen_reader import read_bgen; bgen = read_bgen('merged.bgen', verbose=False)"
       """
        }
        //--variantFilterList snpsToTest.txt.gz --mafFilter 0.01 --hweFilter 1e-06 --callRateFilter 0.95 --machR2Filter 0.4
    //${projectDir}/tools/plink2 --vcf $vcf_file dosage=DS --export bgen-1.2 ref-first --out chr${chr} --extract snpsToTest.txt.gz --chr $chr --maf 0.01 --hwe 1e-06 --geno 0.05 --mac 10  --extract-if-info "R2 > 0.4"
}

/*
 * Convert VCF files to plink and apply variant filtering (per chromosome)
 */
process ConvertVcfToPlink {
    label "medium1"

    input:
        tuple val(chr), path(vcf_file)

    output:
        tuple path("chr*bed"), path("chr*bim"), path("chr*fam"), emit: bfile_per_chr_ch
    
    script:
        if (params.qtls_to_test == ''){
           """
               ${projectDir}/tools/plink2  \
                   --vcf $vcf_file \
                   --make-bed --out chr${chr} \
                   --chr $chr  \
                   --const-fid \
                   --maf 0.01 --hwe 1e-06 --geno 0.05 --mac 10 \
                   --extract-if-info "R2 > 0.4"

           """

        }
        else {
         """
            zcat  ${params.qtls_to_test} | cut -f 2 | gzip > snpsToTest.txt.gz

            ${projectDir}/tools/plink2  \
                --vcf $vcf_file \
                --make-bed --out chr${chr} \
                --chr $chr  \
                --const-fid \
                --extract snpsToTest.txt.gz \
                --maf 0.01 --hwe 1e-06 --geno 0.05 --mac 10 \
                --extract-if-info "R2 > 0.4"

            """
        }

}



/*
 * Combine per chromosome plink genotype files into one
 */
process MergeBgenPerChr {
    label "medium2"

    //echo true
    input:
        path(plink_files)
    output:
        tuple path("merged.bgen"), emit: bfile_ch
    shell:
    '''
        echo !{plink_files}
        exit 1
    '''
}



/*
 * Prepare covariate table: run deconvolution, estimate RNA quality, combine all covariates together
 */
    workflow PREPARE_COVARIATES {
    take:
        exp_type
	    raw_expression_data
        normalized_expression_data
	    signature_matrix_name
	    deconvolution_method
        covariates
	    gene_lengths
        limix_annotation
	    genotype_pcs
	    gte
	    filt_exp_ch

    main:
    //signature matrix path
	signature_matrix = "$projectDir/data/signature_matrices/" + signature_matrix_name  + "_ensg.txt.gz"

    //TODO
	eigenAndIc_ch = MapEigenvectorsAndIcs(filt_exp_ch, params.expression_eigenvectors, params.expression_ics)
	//eigenAndIc_ch = channel.fromPath( '/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/mappedEigenvectorsAndIcs.txt' )
	rnaquality_ch = CalculateRNAQualityScore(normalized_expression_data)
    
    // if no deconvolution is required
    if (deconvolution_method == "NA") {

        CombineCovariatesRNAqual(covariates, Channel.fromPath("NA"), genotype_pcs, gte, rnaquality_ch, eigenAndIc_ch, filt_exp_ch)
        covariates_ch = CombineCovariatesRNAqual.out.covariates_ch
    
    // if lab-based cell proportions are provided: don't run deconvolution
    } else if (deconvolution_method == "lab"){

        cell_counts_ch = Channel.fromPath(params.lab_cell_perc)
        CombineCovariatesRNAqual(covariates, cell_counts_ch, genotype_pcs, gte, rnaquality_ch, eigenAndIc_ch, filt_exp_ch)
        covariates_ch = CombineCovariatesRNAqual.out.covariates_ch

    // the most common case: run deconvolution on TPM-normalized expression, estimate RNA quality and combine all covariates
    } else {
        if (exp_type == "RNAseq" || exp_type == "RNAseq_HGNC") {
            Deconvolution(TPM(raw_expression_data, gene_lengths), signature_matrix, deconvolution_method, exp_type)

        } else {
            Deconvolution(normalized_expression_data, signature_matrix, deconvolution_method, exp_type)
         }
        cell_counts_ch = Deconvolution.out.cellCounts
        CombineCovariatesRNAqual(covariates,cell_counts_ch, genotype_pcs, gte, rnaquality_ch, eigenAndIc_ch, filt_exp_ch)
        covariates_ch = CombineCovariatesRNAqual.out.covariates_ch
    }   
    
    emit:
        covariates_ch

}



