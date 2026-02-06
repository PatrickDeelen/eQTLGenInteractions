#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage() {
  log.info"""
  =======================================================
    Interaction analysis v${workflow.manifest.version}
  =======================================================

  """.stripIndent()
}


