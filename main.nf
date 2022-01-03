#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/ssds
========================================================================================
    Github : https://github.com/nf-core/ssds
    Website: https://nf-co.re/ssds
    Slack  : https://nfcore.slack.com/channels/ssds
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.bwa   = WorkflowMain.getGenomeAttribute(params, 'bwa')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { SSDS } from './workflows/ssds'

//
// WORKFLOW: Run main nf-core/ssds analysis pipeline
//
workflow NFCORE_SSDS {
    SSDS ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_SSDS ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
