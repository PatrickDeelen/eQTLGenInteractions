#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process CreateFiles {

    input:
        val filename
    output:
        path "zScore_*.txt.gz", emit: real
        path "perm/*_perm.txt.gz", emit: permuted

    //[a-z,A-Z,0-9.]
    shell:
        '''

            if [ "!{filename}" = "zScore_ENSG00000100106" ]; then
                echo "c" | gzip > !{filename}A.txt.gz
                echo "d" | gzip > !{filename}A_perm.txt.gz
            fi

            if [ "!{filename}" = "zScore_ENSG00400123123.2" ]; then
                echo "c" | gzip > zScore_ENSG00000100106.txt.gz
                echo "d" | gzip > zScore_ENSG00000100106_perm.txt.gz
            fi


           echo "a" | gzip > !{filename}.txt.gz
           echo "b" | gzip > !{filename}_perm.txt.gz

           mkdir -p ./perm
           mv *_perm.txt.gz ./perm || echo "No results"

        '''

}

process ConcatFiles {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true, failOnError: true

    input :
        path zscoresTest
        path zscorePerm
    output:
        path "result.txt"


    shell:
        '''
            echo !{params.outdir}
            zcat * > result.txt


        '''


}


workflow {


    ch_files = CreateFiles(Channel.of('zScore_ENSG00000100106', 'zScore_ENSG00400123123.2', 'zScore_ENSG00000100104', 'zScore_ENSG00000100106'))

 //    ch_files.real.flatten().unique().view()
//         ch_files.permuted.flatten().unique().view()


    // .map{tuple(it.baseName,it)}.set{ ch_test }


//     concatFiles(  ch_files.real.flatten().unique().collect(),ch_files.permuted.flatten().unique().collect())

        zscoreTest_ch = ch_files.real.flatten().map{it -> [ it.baseName, it ]}.groupTuple().map{ it[1][0]}
        zscorePerm_ch = ch_files.permuted.flatten().map{it -> [ it.baseName, it ]}.groupTuple().map{ it[1][0]}



       //zscoreTest_ch.view()

    ConcatFiles(zscoreTest_ch.collect(),zscorePerm_ch.collect())
}