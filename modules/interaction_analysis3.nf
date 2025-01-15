#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * run QTL mapping per SNP-Gene pair or for all SNPs around the gene depending on command line parameters.
 * Covariates will be included in the model as linear terms, and the $covariate_to_test will also be included as an interaction with genotype
 */
process IeQTLmapping {
    tag "Chunk: $chunk"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true, failOnError: true, pattern: 'beta*'

    input:
       tuple path(tmm_expression), path(covariates), path(limix_annotation), val(chunk), path(qtl_ch)
       path(bgen_ch)
    
    // make the output optional for the case when there are no eQTLs to test and the output is empty. If it's not optional then .collect() in the workflow description will not work
    output:
        path "limix_out/zScore_*.txt.gz", optional: true, emit: zscoresTest
        path "limix_out/perm/zScore_*.txt.gz", optional: true, emit: zscoresPerm
        path "limix_out/feature_metadata_*.txt", optional: true, emit: featureMeta
        path "limix_out/snp_metadata_*.txt", optional: true, emit: snpMeta
        path "limix_out/beta*" , optional: true, emit: betas

    shell:
    '''

    echo !{chunk}

    echo "2"

    outdir=${PWD}/limix_out/
    mkdir -p $outdir

    HOME="./home"
    mkdir -p ${HOME}

    python /groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/singularity_img/Limix_TMP/specialized_lm_interaction_QTL_runner.py \
     --bgen merged \
      -af !{limix_annotation} \
      -cf !{covariates} \
      -pf !{tmm_expression} \
      -od ${outdir} \
      -gr !{chunk} \
      -np !{params.num_perm} \
      -maf 0.01 \
      -hwe 0.0001 \
      --write_permutations \
      --write_zscore \
      -fvf !{qtl_ch} \
      --interaction_term bazinga

    mkdir -p ${outdir}/perm
    mv ${outdir}/*_perm.txt.gz ${outdir}/perm || echo "No results"
      
    '''
}



/*
 * Plot the interaction of rs1981760-NOD2 cis-eQTL with STX3 gene expression (a proxy for neutrophil percentage) - just a sanity check
 */
process PlotSTX3NOD2 {
    label "short"
    publishDir params.outdir, mode: 'copy', overwrite: true, failOnError: true

    input:
        tuple path (expression_path), path (covariate_path), path (bed), path (bim), path (fam), path (exp_PCs)

    output:
        path ("*.pdf")

    shell:
    if (params.expr_pcs == ''){   
        '''
        geno=!{bed}
        plink_base=${geno%.bed}
        echo !{bed}
            !{projectDir}/tools/plink --bfile $plink_base --snp rs1981760 --recode 12 --out snp_geno
        Rscript !{projectDir}/bin/plot_STX3_NOD2.R -e !{expression_path} -c !{covariate_path} -b snp_geno.ped
        '''
    } else {
        '''
        geno=!{bed}
        plink_base=${geno%.bed}
        !{projectDir}/tools/plink --bfile $plink_base --snp rs1981760 --recode 12 --out snp_geno
        Rscript !{projectDir}/bin/plot_STX3_NOD2.R -e !{expression_path} -c !{covariate_path} -b snp_geno.ped -p !{exp_PCs}
        '''
    }
}


/*
 * Convert the interaction analysis output folder to a text file
 */
process ConvertIeQTLsToText {
    debug true

     publishDir params.outdir, mode: 'copy', overwrite: true, failOnError: true

    input:
    path limix_out_files
    
    output:
        path "${params.covariate_to_test}.iqts.txt.gz*"
        //path "limix_out_text/*txt.gz"
        //path "limix_out_text/limix_out_text.md5"


    script:
    """

    python /limix_qtl/Limix_QTL/post_processing/minimal_interaction_postprocess.py \
      -id ./ \
      -od  ./ \
      -sfo \
      --write_compressed

    mv iqtl_results_all.txt.gz ${params.covariate_to_test}.iqts.txt.gz
    md5sum  ${params.covariate_to_test}.iqts.txt.gz > ${params.covariate_to_test}.iqts.txt.gz.md5


    """
}

process CombineResults {
    tag("Merge")

     publishDir params.outdir, mode: 'copy', overwrite: true, failOnError: true

    input:
         path zscoresTest
         path zscorePerm
         path featureMeta
         path snpMeta


    output:
        path "feature_metadata.txt.gz*"
        path "snp_metadata.txt.gz*"
        path "interactionZscoreTest*"
        path "interactionZscorePermutation*"

    script:
    """

        echo -e "feature_id\tchromosome\tstart\tend\tGeneNamebiotype\tn_samples\tn_e_samples"
        tail -n +2 feature_metadata* >> feature_metadata.txt
        gzip feature_metadata.txt
        md5sum feature_metadata.txt > feature_metadata.txt.md5

        echo -e "snp_id\tchromosome\tposition\tassessed_allele\tcall_rate\tmaf\thwe_p"
        tail -n +2 snp_metadata* >> snp_metadata.txt
        gzip snp_metadata.txt
        md5sum snp_metadata.txt > snp_metadata.txt.md5

        java -jar /groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/Datg-tool-1.1/Datg-tool.jar \
            --mode ROW_CONCAT \
            --input ./ \
            --output interactionZscoreTest \
            --filePattern "zScore_([^_]+)\\.txt\\.gz" \
            --datasetName "Interaction z-scores ${params.cohort}" \
            --rowContent "Variant_PermutationRound_eQtlGene" \
            --colContent "Covariates"

        java -jar /groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/Datg-tool-1.1/Datg-tool.jar \
            --mode ROW_CONCAT \
            --input ./ \
            --output interactionZscorePermutation \
            --filePattern "zScore_(.+)_perm\\.txt\\.gz" \
            --datasetName "Permuted interaction z-scores ${params.cohort}" \
            --rowContent "Variant_eQtlGene" \
            --colContent "Covariates"

        echo "Conversion complete"

    """

}

workflow RUN_INTERACTION_QTL_MAPPING {   
    take:
        tmm_expression
        covariates_ch
	    limix_annotation
        covariate_to_test
	    chunk
        qtl_ch
        bgen_ch

    main:

    expr_pcs_ch = params.expr_pcs
            ? Channel.fromPath(params.expr_pcs, checkIfExists:true)
            : Channel.fromPath('EMPTY')
    // if run interaction analysis with covariate * genotype interaction terms, preadjust gene expression for other, linear covariates


       interaction_ch = tmm_expression.combine(covariates_ch).combine(limix_annotation).combine(chunk).combine(qtl_ch)
       IeQTLmapping(interaction_ch, bgen_ch)

       zscoreTest_ch = IeQTLmapping.out.zscoresTest.flatten().unique().collect()
       zscorePerm_ch = IeQTLmapping.out.zscoresPerm.flatten().unique().collect()
       featureMeta_ch = IeQTLmapping.out.featureMeta.flatten().unique().collect()
       snpMeta_ch = IeQTLmapping.out.snpMeta.flatten().unique().collect()


       // zscoreTest_ch.concat(zscorePerm_ch, featureMeta_ch, snpMeta_ch).collect().view()

       CombineResults(zscoreTest_ch, zscorePerm_ch, featureMeta_ch, snpMeta_ch)

    // Plot the interaction of NOD2 cis-eQTL with STX3 (neutrophil proxy) as a sanity check
    //  PlotSTX3NOD2(tmm_expression.combine(covariates_ch).combine(plink_geno).combine(expr_pcs_ch))
}
