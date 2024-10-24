#! /usr/bin/env nextflow

/*
==================================================
                    MetaboXcan
==================================================
 PredictDb Analysis Pipeline.
 #### Homepage / Documentation
 https:://github.com/hakyimlab/metaboxcan
---------------------------------------------------
*/

// Enable DSL 2 syntax
nextflow.enable.dsl = 2


def helpMessage(){
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run main.nf

    Mandatory arguments:
      --gene_models_folder [path]
      --metabolite_models [path]
      --summary_stats_folder [path]
      --covariance [file]
      --trait_name [char]

    Options:
      --gene_models [db]
      --summary_stats [file]

      --zscore_col [column]
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

gene_models = Channel.fromFilePairs(params.gene_models_folder)
//gene_models.view()

metabolite_models = Channel.fromFilePairs(params.metabolite_models_folder)

gwas = Channel
            .fromPath(params.gwas_file)
            .map {gwas_file -> tuple(params.gwas_name,gwas_file)}
//gwas.view()
gene_annot = Channel.fromPath(params.gene_annot)
gene2metabo = Channel.fromPath(params.gene2metabo)
mqtl_db = Channel.fromPath(params.mqtl)
ld_blocks = Channel.fromPath(params.ld_blocks)
metabolite_map = Channel.fromPath(params.metabolite_metadata)

// Import modules
include { spredixcan;smultixcan;smetaboxcan;gwas_database;html_report } from './modules/metaxcan.nf'

// Create workflows
workflow OMICS_PIPELINE {
   take:
      gene_models
      metabolite_models
      gwas
      gene2metabo
      gene_annot
      mqtl_db
      metabolite_map
      ld_blocks

   main:
      // Run Summary prediXcan
      spredixcan(gene_models,gwas.first())

      // Run smultixcan
      smultixcan(spredixcan.out.db.collect(),spredixcan.out.outname.collect(),gwas.first())

      // Run Summary metaboXcan
      smetaboxcan(smultixcan.out,metabolite_models,gwas.first())

      // Convert sumstat to db
      gwas_database(gwas)

      // Generate a report with a locus zoom plot
      html_report(spredixcan.out.outname.collect(), smultixcan.out.multixcan, smetaboxcan.out.smetaboxcan,
                  metabolite_models,gwas_database.out,gene2metabo,gene_annot,mqtl_db,metabolite_map,ld_blocks)
}

workflow {
    main:
      //gwas.view()
      OMICS_PIPELINE(gene_models,metabolite_models,gwas,gene2metabo,gene_annot,mqtl_db,metabolite_map,ld_blocks)
}
