#!/usr/bin/env nextflow

nextflow.enable.dsl = 2 
params.publish_dir_mode = 'copy'
params.samplesheet = false
params.enable_conda = false
params.collect = false
params.skip_flye = false
params.skip_alignments = false
params.flye_mode = '--nano-hq'
params.polish_pilon = false
params.medaka_model = 'r1041_e82_400bps_hac_v4.2.0'
params.skip_ragtag = false
params.out = './results'
/*
 * Print very cool text and parameter info to log. 
 */
log.info """\
==========================================================================================
==========================================================================================
 ▄▄▄       ██▀███   ▄▄▄        ██████   ██████ ▓█████  ███▄ ▄███▓ ▄▄▄▄    ██▓   ▓██   ██▓
▒████▄    ▓██ ▒ ██▒▒████▄    ▒██    ▒ ▒██    ▒ ▓█   ▀ ▓██▒▀█▀ ██▒▓█████▄ ▓██▒    ▒██  ██▒
▒██  ▀█▄  ▓██ ░▄█ ▒▒██  ▀█▄  ░ ▓██▄   ░ ▓██▄   ▒███   ▓██    ▓██░▒██▒ ▄██▒██░     ▒██ ██░
░██▄▄▄▄██ ▒██▀▀█▄  ░██▄▄▄▄██   ▒   ██▒  ▒   ██▒▒▓█  ▄ ▒██    ▒██ ▒██░█▀  ▒██░     ░ ▐██▓░
 ▓█   ▓██▒░██▓ ▒██▒ ▓█   ▓██▒▒██████▒▒▒██████▒▒░▒████▒▒██▒   ░██▒░▓█  ▀█▓░██████▒ ░ ██▒▓░
 ▒▒   ▓▒█░░ ▒▓ ░▒▓░ ▒▒   ▓▒█░▒ ▒▓▒ ▒ ░▒ ▒▓▒ ▒ ░░░ ▒░ ░░ ▒░   ░  ░░▒▓███▀▒░ ▒░▓  ░  ██▒▒▒ 
  ▒   ▒▒ ░  ░▒ ░ ▒░  ▒   ▒▒ ░░ ░▒  ░ ░░ ░▒  ░ ░ ░ ░  ░░  ░      ░▒░▒   ░ ░ ░ ▒  ░▓██ ░▒░ 
  ░   ▒     ░░   ░   ░   ▒   ░  ░  ░  ░  ░  ░     ░   ░      ░    ░    ░   ░ ░   ▒ ▒ ░░  
      ░  ░   ░           ░  ░      ░        ░     ░  ░       ░    ░          ░  ░░ ░     
                                                                       ░         ░ ░     
Niklas Schandry                                                          ░        ░
==========================================================================================
==========================================================================================
  Parameters:
     samplesheet     : ${params.samplesheet}
     collect         : ${params.collect}
     flye_mode       : ${params.flye_mode}
     medaka_model    : ${params.medaka_model}
     polish_pilon    : ${params.polish_pilon}
   Steps skipped:
     skip_flye       : ${params.skip_flye}
     skip_alignments : ${params.skip_alignments}
     skip_ragtag     : ${params.skip_ragtag}
   outdir            : ${params.out}
"""
    .stripIndent(false)

/*
 ===========================================
 * Import processes from modules
 ===========================================
 */

include { COLLECT_READS } from './modules/local/collect_reads/main'
include { ALIGN_TO_BAM as ALIGN } from './modules/align/main'
include { ALIGN_SHORT_TO_BAM as ALIGN_SHORT } from './modules/align/main'
include { BAM_INDEX_STATS_SAMTOOLS as BAM_STATS } from './modules/bam_sort_stat/main'
include { FLYE } from './modules/flye/main'    
include { QUAST } from './modules/quast/main'
include { MEDAKA } from './modules/medaka/main'
include { PILON } from './modules/pilon/main'
include { RAGTAG_SCAFFOLD} from './modules/local/ragtag/main'

/* 
 ===========================================
 ===========================================
 * SUBWORKFLOWS
 ===========================================
 ===========================================
 */

/*
 * COLLECT
 ===========================================
 * Collect reads into single fastq file
 */

workflow COLLECT {
  take: ch_input

  main:
    in_reads = ch_input.map(row -> [row.sample, row.readpath])
    if(params.collect) {
      COLLECT_READS(in_reads)
      in_reads = COLLECT_READS.out.combined_reads
    }
  
  emit:
    in_reads
 }

/*
 * FLYE_ASSEMBLY
 ===========================================
 * Assemble using Flye
 */

workflow FLYE_ASSEMBLY {
  take: in_reads

  main:
    // Asssembly using FLYE
    ch_flye_assembly = Channel.empty()
    FLYE(in_reads, params.flye_mode)
    ch_flye_assembly = FLYE.out.fasta

  emit: ch_flye_assembly
}
/*
 ===========================================
 * MAPPING WORKFLOWS
 ===========================================
 */
/*
 * MAP_TO_REF
 ===========================================
 * map to reference genome using minimap2
 */

workflow MAP_TO_REF {
  take: 
    in_reads
    ch_refs

  main:
    // Map reads to reference
    ch_map_ref_in = in_reads
              .join(ch_refs)
    ALIGN(ch_map_ref_in)
    ch_aln_to_ref = ALIGN.out.alignment
    BAM_STATS(ch_aln_to_ref)

  emit:
    ch_aln_to_ref
}

/*
 * MAP_TO_ASSEMBLY
 ===========================================
 * map to flye assembly using minimap2
 */

workflow MAP_TO_ASSEMBLY {
  take:
    in_reads
    genome_assembly

  main:
    // map reads to assembly
    map_assembly = in_reads
                   .join(genome_assembly) 
    ALIGN(map_assembly)
    aln_to_assembly_bam = ALIGN.out.alignment
    BAM_STATS(aln_to_assembly_bam)
    aln_to_assembly_bai = BAM_STATS.out.bai
    aln_to_assembly_bam_bai = aln_to_assembly_bam.join(aln_to_assembly_bai)

  emit:
    aln_to_assembly_bam
    aln_to_assembly_bai
    aln_to_assembly_bam_bai
}

workflow MAP_SR_TO_ASSEMBLY {
  take:
    in_reads
    genome_assembly

  main:
    // map reads to assembly
    map_assembly = in_reads
                   .join(genome_assembly) 
    ALIGN_SHORT(map_assembly)
    aln_to_assembly_bam = ALIGN_SHORT.out.alignment
    BAM_STATS(ch_aln_to_assembly)
    aln_to_assembly_bai = BAM_STATS.out.bai
    aln_to_assembly_bam_bai = ch_aln_to_assembly.join(aln_to_assembly_bai)

  emit:
    aln_to_assembly_bam
    aln_to_assembly_bai
    aln_to_assembly_bam_bai
}

/*
 * MAP_TO_POLISHED
 ===========================================
 * map to output from polisher using minimap2


workflow MAP_TO_POLISHED {
  take:
    in_reads
    genome_assembly

  main:
    // Remap reads to polished assembly
    map_assembly = in_reads
                   .join(genome_assembly) 
    ch_aln_to_assembly = Channel.empty()
    ALIGN(map_assembly)
    aln_to_assembly_bam = ALIGN.out.alignment
    BAM_STATS(ch_aln_to_assembly)
    aln_to_assembly_bai = BAM_STATS.out.bai
    aln_to_assembly_bam_bai = ch_aln_to_assembly.join(aln_to_assembly_bai)

  emit:
    aln_to_assembly_bam
    aln_to_assembly_bai
    aln_to_assembly_bam_bai
}
*/

/*
 ===========================================
 * QUAST
 ===========================================
 * Run QUAST
 */

workflow RUN_QUAST {
  take: 
    flye_assembly
    ch_input
    aln_to_ref
    aln_to_assembly

  main:
    /* prepare for quast:
     * This makes use of the input channel to obtain the reference and reference annotations
     * See quast module for details
     */
    ch_input_references = ch_input.map(row -> [row.sample, row.ref_fasta, row.ref_gff])
    quast_in = flye_assembly
               .join(ch_input_references)
               .join(aln_to_ref)
               .join(aln_to_assembly)
    /*
     * Run QUAST
     */
    QUAST(quast_in, use_gff = true, use_fasta = true)
}

/*
 ===========================================
 * POLISHING STEPS
 ===========================================
 */

 /*
 * PILON
 ===========================================
 * Run pilon
 */

workflow RUN_PILON {
    take:
      assembly_in
      aln_to_assembly_bam_bai

    main:
      pilon_in = assembly_in.join(aln_to_assembly_bam_bai)
      PILON(pilon_in, "bam")
    
    emit:
      PILON.out.improved_assembly
}

workflow POLISH_PILON {
     take:
       input_channel
       assembly
       ch_aln_to_ref

      main:
          ch_shortreads = input_channel.map { create_shortread_channel(it) }
          MAP_SR_TO_ASSEMBLY(ch_shortreads)
          RUN_PILON(assembly, MAP_SR_TO_ASSEMBLY.out.aln_to_assembly_bam_bai)
          pilon_improved = RUN_PILON.out
          MAP_TO_ASSEMBLY(in_reads, pilon_improved)
          RUN_QUAST(pilon_improved, input_channel, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)
      
      emit:
        pilon_improved
}

 /*
 Accessory function to create input for pilon
 modified from nf-core/rnaseq/subworkflows/local/input_check.nf
 */

def create_shortread_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id       = row.sample
    meta.paired   = row.paired.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def shortreads = []
    if (!file(row.shortread_F).exists()) {
        exit 1, "ERROR: shortread_F fastq file does not exist!\n${row.shortread_F}"
    }
    if (!meta.paired) {
        shortreads = [ meta.id, [ file(row.shortread_F) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: shortread_R fastq file does not exist!\n${row.shortread_R}"
        }
        shortreads = [ meta.id, [ file(row.shortread_F), file(row.shortread_R) ] ]
    }
    return shortreads
}

/*
 * MEDAKA
 ===========================================
 * Run medaka
 */

workflow RUN_MEDAKA {
  take:
     in_reads
     assembly
  
  main:
      medaka_in = in_reads.join(assembly)
      MEDAKA(medaka_in, params.medaka_model)
      medaka_out = MEDAKA.out.assembly
  
  emit:
      medaka_out
}

workflow POLISH_MEDAKA {
    take:
     input_channel
     in_reads
     pilon_improved
     ch_aln_to_ref

    main:
      RUN_MEDAKA(in_reads, pilon_improved)
      medaka_assembly = RUN_MEDAKA.out
      MAP_TO_ASSEMBLY(in_reads, medaka_assembly)
      RUN_QUAST(medaka_assembly, input_channel, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)
    
    emit:
      medaka_assembly
}


/*
 * SCAFFOLDING
 ===========================================
 * Run ragtag scaffold
 */

 workflow RUN_RAGTAG {
  take:
     input_channel
     in_reads
     assembly
     references
     ch_aln_to_ref
  
  main:
      ragtag_in = assembly.join(references)
      RAGTAG_SCAFFOLD(ragtag_in)
      ragtag_scaffold_fasta = RAGTAG_SCAFFOLD.out.corrected_assembly
      ragtag_scaffold_agp = RAGTAG_SCAFFOLD.out.corrected_agp
      MAP_TO_ASSEMBLY(in_reads, ragtag_scaffold_fasta)
      RUN_QUAST(ragtag_scaffold_fasta, input_channel, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)

  emit:
      ragtag_scaffold_fasta
      ragtag_scaffold_agp
}


/*
 ====================================================
 ====================================================
                 MAIN WORKFLOW
 ====================================================
 * Collect fastq files
 * Assemble using flye
      (Enter here: skip_flye)
 * Align to reference 
 * Aling to flye assembly
      (Enter here: skip_flye, skip_alignments)
 * Run quast on flye assembly 
 * Polish with medaka
    * Polish with medaka
    * Align to polished assembly
    * Run quast
  * Polish with pilon (only if polish_pilon true)
    * Polish medaka output with shortreads
    * Align long reads to polished assembly
    * Run quast
 * Scaffold with ragtag
 ====================================================
 ====================================================
 */

workflow {
 /*
 Define channels
 */

 ch_input = Channel.empty()
 ch_refs = Channel.empty()
 ch_aln_to_ref = Channel.empty()
 ch_assembly = Channel.empty()
 ch_assembly_bam = Channel.empty()
 ch_assembly_bam_bai = Channel.empty()
 ch_medaka_in = Channel.empty()
 ch_polished_genome = Channel.empty()
 ch_short_reads = Channel.empty()

  /*
  Check samplesheet
  */

  if(params.samplesheet) {
    ch_input = Channel.fromPath(params.samplesheet) 
                      .splitCsv(header:true) 
    ch_refs = ch_input.map(row -> [row.sample, row.ref_fasta])
                      }
  else {
    exit 1, 'Input samplesheet not specified!'
  }

  /*
  Prepare reads
  */

  COLLECT(ch_input)

  /*
  Prepare assembly
  */

  if(params.skip_flye ) {
    // Sample sheet layout when skipping FLYE
    // sample,readpath,assembly,ref_fasta,ref_gff
    ch_assembly = ch_input.map(row -> [row.sample, row.assembly])
  } else {
    FLYE_ASSEMBLY(COLLECT.out)
    ch_assembly = FLYE_ASSEMBLY.out
  }

  /*
  Prepare alignments
  */

  if(params.skip_alignments) {
    // Sample sheet layout when skipping FLYE and mapping
    // sample,readpath,assembly,ref_fasta,ref_gff,assembly_bam,assembly_bai,ref_bam
    ch_ref_bam = ch_input.map(row -> [row.sample, row.ref_bam]) 
    ch_assembly_bam = ch_input.map(row -> [row.sample, row.assembly_bam]) 
    ch_assembly_bam_bai = ch_input.map(row -> [row.sample, row.assembly_bam, row.assembly_bai]) 
  } else {
    MAP_TO_REF(COLLECT.out, ch_refs)
    ch_ref_bam = MAP_TO_REF.out  

    MAP_TO_ASSEMBLY(COLLECT.out, ch_assembly)
    ch_assembly_bam = MAP_TO_ASSEMBLY.out.aln_to_assembly_bam
    ch_assembly_bam_bai = MAP_TO_ASSEMBLY.out.aln_to_assembly_bam_bai
  }

  /*
  Run QUAST on initial assembly
  */
  RUN_QUAST(ch_assembly, ch_input, ch_ref_bam, ch_assembly_bam)

  /*
  Polishing with medaka
  */

  ch_medaka_in = ch_assembly
  POLISH_MEDAKA(ch_input, COLLECT.out, ch_medaka_in, ch_ref_bam)
  ch_polished_genome = POLISH_MEDAKA.out

  /*
  Polishing with short reads using pulon
  */

  if(params.polish_pilon) {
    POLISH_PILON(ch_input, COLLECT.out, ch_assembly, ch_assembly_bam_bai, ch_ref_bam)
    ch_polished_genome = POLISH_PILON.out.pilon_improved
  } 

  /*
  Scaffolding with ragtag
  */

  if(!params.skip_ragtag) {
    RUN_RAGTAG(ch_input, COLLECT.out, ch_polished_genome, ch_refs, ch_ref_bam)
  }

} 