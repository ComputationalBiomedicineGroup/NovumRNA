#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/novumrna
========================================================================================
    Github : https://github.com/nf-core/novumrna
    Website: https://nf-co.re/novumrna
    Slack  : https://nfcore.slack.com/channels/novumrna
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

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

include { NOVUMRNA } from './workflows/novumrna'

//
// WORKFLOW: Run main nf-core/novumrna analysis pipeline
//
workflow NFCORE_NOVUMRNA {
    NOVUMRNA ()
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
    NFCORE_NOVUMRNA ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
