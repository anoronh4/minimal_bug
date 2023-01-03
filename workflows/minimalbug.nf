/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMinimalbug.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { UCSC_GTFTOGENEPRED             } from '../modules/nf-core/ucsc/gtftogenepred/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST        } from '../modules/nf-core/gatk4/bedtointervallist/main'
include { PICARD_COLLECTRNASEQMETRICS    } from '../modules/nf-core/picard/collectrnaseqmetrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

fasta = "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta"
gtf   = "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.gtf"
bed   = "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.blacklist_intervals.bed"


// Info required for completion email and summary

workflow MINIMALBUG {

    ch_versions = Channel.empty()
     
    bam_ch1 = Channel.of([ [id:'test1',single_end:false ], file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam")])
    bam_ch2 = Channel.of([ [id:'test2',single_end:false ], file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.rna.paired_end.chr6.bam")])
    bam_ch = bam_ch1.mix(bam_ch2).view()

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    GATK4_CREATESEQUENCEDICTIONARY.out.dict.view()
    UCSC_GTFTOGENEPRED(Channel.of([[id:'hg19'],gtf]))
    
    GATK4_BEDTOINTERVALLIST([[id:'hg19'],bed],GATK4_CREATESEQUENCEDICTIONARY.out.dict)

    PICARD_COLLECTRNASEQMETRICS(
        bam_ch,
        UCSC_GTFTOGENEPRED.out.refflat.map{ it -> it[1] },
        fasta,
        GATK4_BEDTOINTERVALLIST.out.interval_list.map{ it -> it[1] }

    )


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
