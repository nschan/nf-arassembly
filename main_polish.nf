#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.publish_dir_mode = 'copy'
params.samplesheet = false
params.enable_conda = false
params.collect = true
params.flye_mode = "--nano-hq"
params.out = './results'
/*
 * Print very cool text and paramter info to log. 
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
==========================================================================================
==========================================================================================
  Parameters:
   samplesheet   : ${params.samplesheet}
   collect       : ${params.collect}
   flye_mode     : ${params.flye_mode}
   outdir        : ${params.out}
"""
    .stripIndent(false)
/*
 * Import processes from modules
 */
include { COLLECT_READS } from './modules/local/collect_reads/main'
include { ALIGN_TO_BAM as ALIGN } from './modules/align/main'
include { BAM_INDEX_STATS_SAMTOOLS as BAM_STATS } from './modules/bam_sort_stat/main'
include { FLYE } from './modules/flye/main'    
include { QUAST } from './modules/quast/main'
inlcude { MEDAKA } from './modules/medaka/main'
inlcude { PILON } from './modules/medaka/main'
/* 
 * SUBWORKFLOWS
 */

/*
 * SUBWORKFLOW:
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
 * SUBWORKFLOW:
 * assembly via Flye
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
 * SUBWORKFLOW:
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
 * SUBWORKFLOW:
 * map to flye assembly using minimap2
 */

workflow MAP_TO_ASSEMBLY {
  take:
    in_reads
    genome_assembly

  main:
    // Remap reads to flye assembly
    map_assembly = in_reads
                   .join(genome_assembly) 
    ch_aln_to_assembly = Channel.empty()
    ALIGN(map_assembly)
    ch_aln_to_assembly = ALIGN.out.alignment
    BAM_STATS(ch_aln_to_assembly)
    aln_to_assembly_bai = BAM_STATS.out.bai
    aln_to_assembly_bam_bai = ch_aln_to_assembly.join(aln_to_assembly_bai)

  emit:
    ch_aln_to_assembly
    aln_to_assembly_bam_bai
}

workflow MAP_TO_POLISHED {
  take:
    in_reads
    genome_assembly

  main:
    // Remap reads to flye assembly
    map_assembly = in_reads
                   .join(genome_assembly) 
    ch_aln_to_assembly = Channel.empty()
    ALIGN(map_assembly)
    ch_aln_to_assembly = ALIGN.out.alignment
    BAM_STATS(ch_aln_to_assembly)
    aln_to_assembly_bai = BAM_STATS.out.bai
    aln_to_assembly_bam_bai = ch_aln_to_assembly.join(aln_to_assembly_bai)

  emit:
    ch_aln_to_assembly
    aln_to_assembly_bam_bai
}

/*
 * SUBWORKFLOW:
 * Run quast
 */

workflow RUN_QUAST {
  take: 
    flye_assembly
    ch_input
    aln_to_ref
    aln_to_assembly

  main:
    /* prepare for quast
     * This makes use of the input channel to obtain the reference and reference annotations
     * See quast module for details
     */
    ch_input_references = ch_input.map(row -> [row.sample, row.ref_fasta, row.ref_gff])
    quast_in = flye_assembly
               .join(ch_input_references)
               .join(aln_to_ref)
               .join(aln_to_assembly)
    QUAST(quast_in, use_gff = true, use_fasta = true)
}

/*
 * ADD POLISHING STEP
 */
workflow RUN_PILON {
    take:
      flye_assembly
      aln_to_assembly_bam_bai

    main:
      pilon_in = flye_assembly.join(aln_to_assembly_bam_bai)
      PILON(pilon_in, "bam")
    
    emit:
      improved_assembly
}

workflow POLISH_PILON {
     take:
       input_channel
       in_reads
       flye_assembly
       ch_aln_to_ref

      main:
          RUN_PILON(flye_assembly, aln_to_assembly_bam_bai)
          pilon_improved = RUN_PILON.out.improved_assembly
          MAP_TO_POLISHED(in_reads, RUN_PILON.out.improved_assembly)
          RUN_QUAST(RUN_PILON.out, input_channel, ch_aln_to_ref, MAP_TO_POLISHED.out)
      
      emit:
        pilon_improved
    }

workflow RUN_MEDAKA {
  take:
     in_reads
     improved_assembly
  
  main:
      medaka_in = in_reads.join(improved_assembly)
      MEDAKA(medaka_in)
  
  emit: 
     assembly
}

workflow POLISH_MEDAKA {
    take:
     input_channel
     in_reads
     pilon_improved
     ch_aln_to_ref

    main:
      RUN_MEDAKA(in_reads, pilon_improved)
      MAP_TO_POLISHED(in_reads, RUN_MEDAKA.out.assembly)
      RUN_QUAST(RUN_MEDAKA.out, input_channel, ch_aln_to_ref, MAP_TO_POLISHED.out)
}

/*
 * MAIN WORKFLOW
 * Run Collect
 * Run flye
 * Run minimap2 to align to reference and flye assembly
 * Run quast
 */

workflow {
    // Sample sheet layout:
    // sample, readpath, ref_fasta, ref_gff
    ch_input = Channel.fromPath(params.samplesheet) 
                      .splitCsv(header:true)
    ch_refs = ch_input.map(row -> [row.sample, row.ref_fasta])

    COLLECT(ch_input)

    FLYE_ASSEMBLY(COLLECT.out)

    MAP_TO_REF(COLLECT.out, ch_refs)

    MAP_TO_ASSEMBLY(COLLECT.out, FLYE_ASSEMBLY.out)

    RUN_QUAST(FLYE_ASSEMBLY.out, ch_input, MAP_TO_REF.out, MAP_TO_ASSEMBLY.out.ch_aln_to_assembly)

    POLISH_PILON(ch_input, COLLECT.out, FLYE_ASSEMBLY.out, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam_bai, MAP_TO_REF.out)

    POLISH_MEDAKA(ch_input, COLLECT.out, POLISH_PILON.pilon_improved, MAP_TO_REF.out)

    


}