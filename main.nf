#!/usr/bin/env nextflow

nextflow.enable.dsl = 2 
params.publish_dir_mode = 'copy'
params.samplesheet = false
params.enable_conda = false
params.collect = false
params.skip_flye = false
params.skip_alignments = false
params.flye_mode = '--nano-hq'
params.flye_args = ''
params.polish_pilon = false
params.medaka_model = 'r1041_e82_400bps_hac_v4.2.0'
params.scaffold_ragtag = false
params.scaffold_links = false
params.scaffold_slr = false
params.scaffold_longstitch = false
params.lift_annotations = true
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
                                                                          ░        ░
------------------------------------------------------------------------------------------
Niklas Schandry      niklas@bio.lmu.de      https://gitlab.lrz.de/beckerlab/nf-arassembly                                          
------------------------------------------------------------------------------------------
  Results directory  : ${params.out}

  Parameters:
     samplesheet     : ${params.samplesheet}
     collect         : ${params.collect}
     flye_mode       : ${params.flye_mode}
     medaka_model    : ${params.medaka_model}
     polish_pilon    : ${params.polish_pilon}

    Steps skipped
     skip_flye       : ${params.skip_flye}
     skip_alignments : ${params.skip_alignments}

    Scaffolding Tools
      ragtag         : ${params.scaffold_ragtag}
      LINKS          : ${params.scaffold_links}
      SLR            : ${params.scaffold_slr}
      longstitch     : ${params.scaffold_longstitch}

    Annotation lift  : ${params.lift_annotations}

==========================================================================================
==========================================================================================
"""
    .stripIndent(false)

/*
 ===========================================
 * Import processes from modules
 ===========================================
 */

// Preprocessing
include { COLLECT_READS } from './modules/local/collect_reads/main'

// Read statistics
include { NANOQ } from './modules/nanoq/main'

// Alignment
include { ALIGN_TO_BAM as ALIGN } from './modules/align/main'
include { ALIGN_SHORT_TO_BAM as ALIGN_SHORT } from './modules/align/main'
include { BAM_INDEX_STATS_SAMTOOLS as BAM_STATS } from './modules/bam_sort_stat/main'

// Assembly 
include { FLYE } from './modules/flye/main'    

// Polishing
include { MEDAKA } from './modules/medaka/main'
include { PILON } from './modules/pilon/main'

// Scaffolding
include { RAGTAG_SCAFFOLD } from './modules/local/ragtag/main'
include { LINKS } from './modules/local/links/main'
include { SLR } from './modules/local/slr/main'
include { LONGSTITCH } from './modules/local/longstitch/main'

// Annotation
include { LIFTOFF } from './modules/local/liftoff/main'

// Quality control
include { QUAST } from './modules/quast/main'
include { BUSCO } from './modules/busco/main'

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

 workflow RUN_NANOQ {
  take: in_reads

  main:
  
  NANOQ(in_reads)
  report = NANOQ.out.report
  stats = NANOQ.out.stats

  emit:
    report
    stats
 }

 /*
 * FLYE_ASSEMBLY
 ===========================================
 * Assemble using Flye
 */

workflow FLYE_ASSEMBLY {
  take: in_reads
        ch_input

  main:
    // Asssembly using FLYE
    ch_flye_assembly = Channel.empty()
    FLYE(in_reads, params.flye_mode)
    ch_flye_assembly = FLYE.out.fasta
    if(params.lift_annotations) {
    RUN_LIFTOFF(FLYE.out.fasta, ch_input)
  }
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
 * map to assembly using minimap2
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

workflow MAP_SR {
  take:
    in_reads
    genome_assembly

  main:
    // map reads to assembly
    map_assembly = in_reads
                   .join(genome_assembly) 
    ALIGN_SHORT(map_assembly)
    aln_to_assembly_bam = ALIGN_SHORT.out.alignment
    BAM_STATS(aln_to_assembly_bam)
    aln_to_assembly_bai = BAM_STATS.out.bai
    aln_to_assembly_bam_bai = aln_to_assembly_bam.join(aln_to_assembly_bai)

  emit:
    aln_to_assembly_bam
    aln_to_assembly_bai
    aln_to_assembly_bam_bai
}

/*
 ===========================================
 * POLISHING STEPS
 ===========================================
 */

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
     ch_input
     in_reads
     pilon_improved
     ch_aln_to_ref

    main:
      RUN_MEDAKA(in_reads, pilon_improved)
      medaka_assembly = RUN_MEDAKA.out
      MAP_TO_ASSEMBLY(in_reads, medaka_assembly)
      RUN_QUAST(medaka_assembly, ch_input, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)
      RUN_BUSCO(medaka_assembly)
      if(params.lift_annotations) {
        RUN_LIFTOFF(RUN_MEDAKA.out, ch_input)
      }
    emit:
      medaka_assembly
}

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
       ch_input
       in_reads
       assembly
       ch_aln_to_ref

      main:
          ch_shortreads = ch_input.map { create_shortread_channel(it) }
          MAP_SR(ch_shortreads, assembly)
          RUN_PILON(assembly, MAP_SR.out.aln_to_assembly_bam_bai)
          pilon_improved = RUN_PILON.out
          MAP_TO_ASSEMBLY(in_reads, pilon_improved)
          RUN_QUAST(pilon_improved, ch_input, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)
          RUN_BUSCO(pilon_improved)
          if(params.lift_annotations) {
             RUN_LIFTOFF(RUN_PILON.out, ch_input)
          }
      
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
        shortreads = [ meta.id, meta.paired, [ file(row.shortread_F) ] ]
    } else {
        if (!file(row.shortread_R).exists()) {
            exit 1, "ERROR: shortread_R fastq file does not exist!\n${row.shortread_R}"
        }
        shortreads = [ meta.id, meta.paired, [ file(row.shortread_F), file(row.shortread_R) ] ]
    }
    return shortreads
}

 /*
 * SCAFFOLDING
 ===========================================
 * Run ragtag scaffold
 * Run LINKS
 * Run SLR
 */

 workflow RUN_RAGTAG {
  take:
     ch_input
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
      RUN_QUAST(ragtag_scaffold_fasta, ch_input, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)
      RUN_BUSCO(ragtag_scaffold_fasta)
      if(params.lift_annotations) {
             RUN_LIFTOFF(RAGTAG_SCAFFOLD.out.corrected_assembly, ch_input)
      }

  emit:
      ragtag_scaffold_fasta
      ragtag_scaffold_agp
}

 workflow RUN_LINKS {
  take:
     ch_input
     in_reads
     assembly
     references
     ch_aln_to_ref
  
  main:
      links_in = assembly.join(in_reads)
      LINKS(links_in)
      scaffolds = LINKS.out.scaffolds
      MAP_TO_ASSEMBLY(in_reads, scaffolds)
      RUN_QUAST(scaffolds, ch_input, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)
      RUN_BUSCO(scaffolds)
      if(params.lift_annotations) {
             RUN_LIFTOFF(LINKS.out.scaffolds, ch_input)
      }

  emit:
     scaffolds
}

workflow RUN_SLR {
  take:
     ch_input
     in_reads
     assembly
     references
     ch_aln_to_ref
  
  main:
      slr_in = assembly.join(in_reads)
      SLR(slr_in)
      scaffolds = SLR.out.scaffolds
      MAP_TO_ASSEMBLY(in_reads, scaffolds)
      RUN_QUAST(scaffolds, ch_input, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)
      RUN_BUSCO(scaffolds)
      if(params.lift_annotations) {
        RUN_LIFTOFF(SLR.out.scaffolds, ch_input)
      }

  emit:
     scaffolds
}
workflow RUN_LONGSTITCH {
  take:
     ch_input
     in_reads
     assembly
     references
     ch_aln_to_ref
  
  main:
      longstitch_in = assembly.join(in_reads)
      LONGSTITCH(longstitch_in)
      scaffolds = LONGSTITCH.out.scaffolds
      MAP_TO_ASSEMBLY(in_reads, scaffolds)
      RUN_QUAST(scaffolds, ch_input, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)
      RUN_BUSCO(scaffolds)
      if(params.lift_annotations) {
        RUN_LIFTOFF(LONGSTITCH.out.scaffolds, ch_input)
      }

  emit:
     scaffolds
}

/*
 ===========================================
 * QUALITY CONTROL STEPS
 ===========================================
 */

 /*
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
 * BUSCO
 ===========================================
 * Run BUSCO standalone 
 * BUSCO in quast is quite old.
 */

workflow RUN_BUSCO {
  take: 
    assembly

  main:
    busco_in = assembly
    /*
     * Run BUSCO
     */
    BUSCO(busco_in ,"brassicales_odb10", "/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/software/busco_db")
}

workflow RUN_LIFTOFF {
  take:
    assembly
    inputs
  
  main:
    liftoff_in = assembly.join(inputs.map(row -> [row.sample, row.ref_fasta, row.ref_gff]))
    LIFTOFF(liftoff_in)
    lifted_annotations = LIFTOFF.out

  emit:
    lifted_annotations
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
 * Scaffold with ragtag, LINKS, SLR or LONGSTITCH
 ====================================================
 ====================================================
 */

workflow ARASEMMBLY {
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

  NANOQ(COLLECT.out)
  /*
  Prepare assembly
  */

  if(params.skip_flye ) {
    // Sample sheet layout when skipping FLYE
    // sample,readpath,assembly,ref_fasta,ref_gff
    ch_assembly = ch_input.map(row -> [row.sample, row.assembly])
  } else {
    FLYE_ASSEMBLY(COLLECT.out,ch_input)
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
  RUN_BUSCO(ch_assembly)

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
    POLISH_PILON(ch_input, COLLECT.out, ch_polished_genome, ch_ref_bam)
    ch_polished_genome = POLISH_PILON.out.pilon_improved
  } 

  /*
  Scaffolding
  */

  if(params.scaffold_ragtag) {
    RUN_RAGTAG(ch_input, COLLECT.out, ch_polished_genome, ch_refs, ch_ref_bam)
  }

  if(params.scaffold_links) {
    RUN_LINKS(ch_input, COLLECT.out, ch_polished_genome, ch_refs, ch_ref_bam)
  }

  if(params.scaffold_slr) {
    RUN_SLR(ch_input, COLLECT.out, ch_polished_genome, ch_refs, ch_ref_bam)
  }
  
  if(params.scaffold_longstitch) {
    RUN_LONGSTITCH(ch_input, COLLECT.out, ch_polished_genome, ch_refs, ch_ref_bam)
  }
  
} 

 /*
 ====================================================
 ====================================================
                 RUN PIPELINE
 ====================================================
 ====================================================
 */
workflow {
  ARASEMMBLY()
}