#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.publish_dir_mode = 'copy'
params.samplesheet = false
params.enable_conda = false
params.collect = true
params.flye_mode = "--nano-hq"
params.out = './results'

include { COLLECT_READS } from './modules/local/collect_reads/main'
include { ALIGN_TO_BAM as ALIGN } from './modules/align/main'
include { BAM_INDEX_STATS_SAMTOOLS as BAM_STATS } from './modules/bam_sort_stat/main'
include { FLYE } from './modules/flye/main'    
include { QUAST as RUN_QUAST } from './modules/quast/main'

workflow COLLECT {
  take: ch_input

  main:
    in_reads = ch_input.map(row -> [row.sample, row.readpath])
    if(params.collect) {
      // COLLECT
      // Collect reads from fastq.gz into single fastq
  
      // view for debug
      // in_reads.view()
      COLLECT_READS(in_reads)
      in_reads = COLLECT_READS.out.combined_reads
    } else {
      in_reads = in_reads.map{sample, readpath -> [sample, combined_reads]}   
    }

  emit:
    in_reads
 }

workflow FLYE_ASSEMBLY {
  take: in_reads

  main:
    // Asssembly using FLYE
    ch_flye_assembly = Channel.empty()
    FLYE(in_reads, params.flye_mode)
    ch_flye_assembly = FLYE.out.fasta

  emit: ch_flye_assembly
}

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

workflow MAP_TO_ASSEMBLY {
  take:
    in_reads
    flye_assembly

  main:
    // Remap reads to flye assembly
    map_assembly = in_reads
                   .join(flye_assembly) 
    ch_aln_to_assembly = Channel.empty()
    ALIGN(map_assembly)
    ch_aln_to_assembly = ALIGN.out.alignment
    BAM_STATS(ch_aln_to_assembly)

  emit:
    ch_aln_to_assembly
}

workflow QUAST {
  take: 
    flye_assembly
    ch_input
    aln_to_ref
    aln_to_assembly

  main:
    // prepare for quast
    ch_input_references = ch_input.map(row -> [row.sample, row.ref_fasta, row.ref_gff])
    quast_in = flye_assembly
               .join(ch_input_references)
               .join(aln_to_ref)
               .join(aln_to_assembly)
    RUN_QUAST(quast_in, use_gff = true, use_fasta = true)
}

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
    
    QUAST(FLYE_ASSEMBLY.out, ch_input, MAP_TO_REF.out, MAP_TO_ASSEMBLY.out)
}