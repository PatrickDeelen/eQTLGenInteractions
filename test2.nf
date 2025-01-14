#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process createFiles {

    input:
        val filename
    output:
        path "zScore_[a-z,A-Z,0-9]*.txt.gz", emit: real
        path "*_perm.txt.gz", emit: permuted

    //[a-z,A-Z,0-9.]
    shell:
        '''

            if [ "!{filename}" = "zScore_ENSG00000100106" ]; then
                echo "c" | gzip > !{filename}A.txt.gz
                echo "d" | gzip > !{filename}A_perm.txt.gz
            fi

           echo "a" | gzip > !{filename}.txt.gz
           echo "b" | gzip > !{filename}_perm.txt.gz
        '''

}

process concatFiles {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true, failOnError: true

    input :
        path '*'
    output:
        "result.txt"


    shell:
        '''
            echo !{params.outdir}
            zcat * > result.txt


        '''


}


workflow {


    ch_files = createFiles(Channel.of('zScore_ENSG00000100106', 'zScore_ENSG00400123123.2', 'zScore_ENSG00000100104', 'zScore_ENSG00000100106'))

    ch_files.real.flatten().unique().view()

    // .map{tuple(it.baseName,it)}.set{ ch_test }

    // ch_test.view()

    // concatFiles(ch_test)


}