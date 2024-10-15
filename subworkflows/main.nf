/*
 ===========================================
 * Import processes from modules
 ===========================================
 */

// Preprocessing
include { COLLECT_READS } from '../modules/local/collect_reads/main'
include { PORECHOP } from '../modules/porechop/main'
include { LIMA } from '../modules/lima/main'
include { SAMTOOLS_FASTQ as TO_FASTQ } from '../modules/samtools/fastq/main'
// Read statistics
include { NANOQ } from '../modules/nanoq/main'

// Jellyfish
include { COUNT } from '../modules/jellyfish/main'
include { DUMP } from '../modules/jellyfish/main'
include { HISTO } from '../modules/jellyfish/main'
include { STATS } from '../modules/jellyfish/main'
include { GENOMESCOPE } from '../modules/genomescope/main'

// Alignment
include { ALIGN_TO_BAM as ALIGN } from '../modules/align/main'
include { ALIGN_SHORT_TO_BAM as ALIGN_SHORT } from '../modules/align/main'
include { BAM_INDEX_STATS_SAMTOOLS as BAM_STATS } from '../modules/bam_sort_stat/main'

// Assembly 
include { FLYE } from '../modules/flye/main'    
include { HIFIASM; HIFIASM_UL } from '../modules/hifiasm/main'

// Polishing
include { MEDAKA } from '../modules/medaka/main'
include { PILON } from '../modules/pilon/main'

// Scaffolding
include { RAGTAG_SCAFFOLD } from '../modules/local/ragtag/main'
include { LINKS } from '../modules/local/links/main'
include { LONGSTITCH } from '../modules/local/longstitch/main'

// Annotation
include { LIFTOFF } from '../modules/local/liftoff/main'

// Quality control
include { QUAST } from '../modules/quast/main'
include { BUSCO } from '../modules/busco/main'

// YAK
include { KMER_ASSEMBLY } from '../modules/yak/main'
include { KMER_LONGREADS as KMER_ONT } from '../modules/yak/main'
include { KMER_LONGREADS as KMER_HIFI } from '../modules/yak/main'
include { KMER_SHORTREADS } from '../modules/yak/main'
include { KMER_HISTOGRAM as KMER_ONT_HIST } from '../modules/yak/main'
include { KMER_HISTOGRAM as KMER_HIFI_HIST} from '../modules/yak/main'
include { KMER_HISTOGRAM as KMER_SR_HIST} from '../modules/yak/main'
include { KMER_HISTOGRAM as KMER_ASSEMBLY_HIST} from '../modules/yak/main'
include { READ_QV as KMER_ONT_QV } from '../modules/yak/main'
include { READ_QV as KMER_HIFI_QV } from '../modules/yak/main'
include { ASSEMBLY_KQV as KMER_ASSEMBLY_QV } from '../modules/yak/main'

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
  
    ch_input
      .map { row -> [row.sample, row.ontreads] }
      .set { in_reads }

    if(params.collect) {
      COLLECT_READS(in_reads)

      COLLECT_READS
        .out
        .combined_reads
        .set { in_reads }
    }

  emit:
    in_reads
 }
 workflow CHOP {

  take: in_reads

  main:
  
  if(params.porechop) {
    PORECHOP(in_reads)
    PORECHOP
      .out
      .reads
      .set { chopped_reads }
  } else {
    in_reads
      .set { chopped_reads }
  }
  
  
  emit:
    chopped_reads
 }

 workflow RUN_NANOQ {
  take: in_reads

  main:
  
  NANOQ(in_reads)

  NANOQ
    .out
    .report
    .set { report }

  NANOQ
    .out
    .stats
    .set { stats }

  NANOQ
    .out
    .stats
    .set { median_length }

  emit:
    report
    stats
    median_length
 }
 
 /*
 nanopore specific read preparation workfflow
 */

 workflow PREPARE_ONT {
  take: inputs

  main:
    COLLECT(inputs)

    CHOP(COLLECT.out)

    trimmed = CHOP.out

    NANOQ(trimmed)

    trimmed_med_len = NANOQ.out.median_length

  emit:
      trimmed
      trimmed_med_len

 }

 workflow JELLYFISH {
    take:
        samples // id, fasta
        nanoq_out
    
    main: 
        COUNT(samples)
        COUNT
          .out
          .set { kmers }
          
        if(params.dump) {
          DUMP(kmers)
        }

        HISTO(kmers)
        if(!params.read_length == null) {
            HISTO
            .out
            .map { it -> [it[0], it[1], params.kmer_length, params.read_length] }
            .set { genomescope_in }
        } 
        if(params.read_length == null) {
            HISTO
            .out
            .map { it -> [it[0], it[1], params.kmer_length] }
            .join( nanoq_out )
            .set { genomescope_in }
        }
        GENOMESCOPE(genomescope_in)

        STATS(kmers)
        GENOMESCOPE.out.estimated_hap_len
          .set{ hap_len }

    emit: hap_len
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
    in_reads
      .join(ch_refs)
      .set { ch_map_ref_in }

    ALIGN(ch_map_ref_in)

    ALIGN
      .out
      .alignment
      .set { ch_aln_to_ref }

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
    in_reads
      .join(genome_assembly)
      .set { map_assembly }

    ALIGN(map_assembly)

    ALIGN
      .out
      .alignment
      .set { aln_to_assembly_bam }

    BAM_STATS(aln_to_assembly_bam)

    BAM_STATS
      .out
      .bai
      .set { aln_to_assembly_bai }

    aln_to_assembly_bam
      .join(aln_to_assembly_bai)
      .set { aln_to_assembly_bam_bai }

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
    in_reads
      .join(genome_assembly)
      .set { map_assembly }

    ALIGN_SHORT(map_assembly)

    ALIGN_SHORT
      .out
      .alignment
      .set { aln_to_assembly_bam }

    BAM_STATS(aln_to_assembly_bam)

    BAM_STATS
      .out
      .bai
      .set { aln_to_assembly_bai }

    aln_to_assembly_bam.
      join(aln_to_assembly_bai)
      .set { aln_to_assembly_bam_bai }

  emit:
    aln_to_assembly_bam
    aln_to_assembly_bai
    aln_to_assembly_bam_bai
}

/*
 ===========================================
 ASSEMBLY
 ===========================================
*/

workflow ASSEMBLE_ONT {
  take: 
    ont_reads
    ch_input
    genomescope_out
    shortread_kmers
    
  main:
    Channel.empty().set { ch_refs }

    if (params.use_ref) {
      ch_input
        .map { row -> [row.sample, row.ref_fasta] }
        .set { ch_refs }
    }

    if(params.skip_flye ) {
      // Sample sheet layout when skipping FLYE
      // sample,ontreads,assembly,ref_fasta,ref_gff
      ch_input
        .map { row -> [row.sample, row.assembly] }
        .set { ch_assembly }
    } else {
      // Run flye
      if(!params.genome_size == null) {
        ont_reads
          .map { it -> [it[0], it[1], params.genome_size] }
          .set { flye_in }
      }
      if(params.genome_size == null) {
        ont_reads
          .join(genomescope_out)
          .set { flye_in }
      } 
      FLYE(flye_in, params.flye_mode)

      FLYE
        .out
        .fasta
        .set { ch_assembly }
    }

    /*
    Prepare alignments
    */

    if(params.skip_alignments) {
      // Sample sheet layout when skipping FLYE and mapping
      // sample,ontreads,assembly,ref_fasta,ref_gff,assembly_bam,assembly_bai,ref_bam
      ch_input
        .map { row -> [row.sample, row.ref_bam] }
        .set { ch_ref_bam }

      ch_input
        .map { row -> [row.sample, row.assembly_bam] }
        .set { ch_assembly_bam }

      ch_input
        .map { row -> [row.sample, row.assembly_bam, row.assembly_bai] }
        .set { ch_assembly_bam_bai } 

    } else {
      Channel.empty().set { ch_ref_bam }

      if(params.use_ref) {
        MAP_TO_REF(ont_reads, ch_refs)

        MAP_TO_REF
          .out
          .set { ch_ref_bam }
      }

      MAP_TO_ASSEMBLY(ont_reads, ch_assembly)
      MAP_TO_ASSEMBLY
        .out
        .aln_to_assembly_bam
        .set { ch_assembly_bam }

      MAP_TO_ASSEMBLY
        .out
        .aln_to_assembly_bam_bai
        .set { ch_assembly_bam_bai }
    }
    if(params.lift_annotations) {
      RUN_LIFTOFF(FLYE.out.fasta, ch_input)
    }
    /*
    QC on initial assembly
    */
    YAK_QC(ch_assembly, shortread_kmers)

    RUN_QUAST(ch_assembly, ch_input, ch_ref_bam, ch_assembly_bam)

    RUN_BUSCO(ch_assembly)

  emit: 
    ch_assembly
    ch_ref_bam
}

workflow ASSEMBLE_HIFI {
  take: 
    hifi_reads // normal mode: meta, hifireads; UL mode: meta, hifireads, ontreads
    ch_input
    shortread_kmers
        
  main:
    Channel.empty().set { ch_refs }

    if (params.use_ref) {
      ch_input
        .map { row -> [row.sample, row.ref_fasta] }
        .set { ch_refs }
    }

    if(params.hifi_ont) {
      HIFIASM_UL(hifi_reads, params.hifi_args)

      HIFIASM_UL
        .out
        .primary_contigs_fasta
        .set { ch_assembly }
    } else {
      HIFIASM(hifi_reads, params.hifi_args)

      HIFIASM
        .out
        .primary_contigs_fasta
        .set { ch_assembly }
    }

    /*
    Prepare alignments
    */

    if(params.skip_alignments) {
      // Sample sheet layout when skipping FLYE and mapping
      // sample,ontreads,assembly,ref_fasta,ref_gff,assembly_bam,assembly_bai,ref_bam
      ch_input
        .map { row -> [row.sample, row.ref_bam] }
        .set { ch_ref_bam }

      ch_input
        .map { row -> [row.sample, row.assembly_bam] }
        .set { ch_assembly_bam }

      ch_input
        .map { row -> [row.sample, row.assembly_bam, row.assembly_bai] }
        .set { ch_assembly_bam_bai } 

    } else {
      Channel.empty().set { ch_ref_bam }

      if(params.hifi_ont) {
        hifi_reads
          .map { it -> [it[0], it[1]]}
          .set { ch_reads }
      } else {
        hifi_reads
          .set { ch_reads }
      }

      if(params.use_ref) {
        MAP_TO_REF(ch_reads, ch_refs)

        MAP_TO_REF
          .out
          .set { ch_ref_bam }
      }

      MAP_TO_ASSEMBLY(ch_reads, ch_assembly)
      MAP_TO_ASSEMBLY
        .out
        .aln_to_assembly_bam
        .set { ch_assembly_bam }

      MAP_TO_ASSEMBLY
        .out
        .aln_to_assembly_bam_bai
        .set { ch_assembly_bam_bai }
    }
    
    if(params.lift_annotations) {
      RUN_LIFTOFF(ch_assembly, ch_input)
    }

  /*
  QC on initial assembly
  */ 
    YAK_QC(ch_assembly, shortread_kmers)

    RUN_QUAST(ch_assembly, ch_input, ch_ref_bam, ch_assembly_bam)

    RUN_BUSCO(ch_assembly)

  emit:
    ch_assembly
}

/*
 ===========================================
 * POLISHING
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
    in_reads
      .join(assembly)
      .set { medaka_in }

    MEDAKA(medaka_in, params.medaka_model)

    MEDAKA
      .out
      .assembly
      .set { medaka_out }
  
  emit:
     medaka_out
}

workflow POLISH_MEDAKA {
    take:
      ch_input
      in_reads
      assembly
      ch_aln_to_ref
      shortread_kmers

    main:
      RUN_MEDAKA(in_reads, assembly)
      
      MAP_TO_ASSEMBLY(in_reads, RUN_MEDAKA.out)

      RUN_QUAST(RUN_MEDAKA.out, ch_input, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)

      RUN_BUSCO(RUN_MEDAKA.out)

      YAK_QC(RUN_MEDAKA.out, shortread_kmers)

      if(params.lift_annotations) {
        RUN_LIFTOFF(RUN_MEDAKA.out, ch_input)
      }
    emit:
      RUN_MEDAKA.out
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
      assembly_in
        .join(aln_to_assembly_bam_bai)
        .set { pilon_in }

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
    shortread_kmers
  
  main:
    ch_input
        .map { create_shortread_channel(it) }
        .set { ch_shortreads }
    MAP_SR(ch_shortreads, assembly)

    RUN_PILON(assembly, MAP_SR.out.aln_to_assembly_bam_bai)

    RUN_PILON
       .out
       .set { pilon_improved }

    MAP_TO_ASSEMBLY(in_reads, pilon_improved)

    RUN_QUAST(pilon_improved, ch_input, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)
   
    RUN_BUSCO(pilon_improved)

    YAK_QC(pilon_improved, shortread_kmers)


    if(params.lift_annotations) {
        RUN_LIFTOFF(RUN_PILON.out, ch_input)
    }
  
  emit:
    pilon_improved
}

/*
 ===========================================
 * SCAFFOLDING
 ===========================================
 */

/* Workflows for scaffolding tools
 ===========================================
 * Run ragtag scaffold
 * Run LINKS
 * Run SLR
 * Run LONGSTITCH
 */

workflow RUN_RAGTAG {
  take:
    inputs
    in_reads
    assembly
    references
    ch_aln_to_ref
    shortread_kmers
  
  main:
    assembly
      .join(references)
      .set { ragtag_in }

    RAGTAG_SCAFFOLD(ragtag_in)

    RAGTAG_SCAFFOLD
      .out
      .corrected_assembly
      .set { ragtag_scaffold_fasta }

    RAGTAG_SCAFFOLD
      .out
      .corrected_agp
      .set { ragtag_scaffold_agp }

    MAP_TO_ASSEMBLY(in_reads, ragtag_scaffold_fasta)

    RUN_QUAST(ragtag_scaffold_fasta, inputs, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)

    RUN_BUSCO(ragtag_scaffold_fasta)

    YAK_QC(ragtag_scaffold_fasta, shortread_kmers)

    if(params.lift_annotations) {
      RUN_LIFTOFF(RAGTAG_SCAFFOLD.out.corrected_assembly, inputs)
    }

  emit:
      ragtag_scaffold_fasta
      ragtag_scaffold_agp
}

workflow RUN_LINKS {
  take:
    inputs
    in_reads
    assembly
    references
    ch_aln_to_ref
    shortread_kmers
  
  main:
    assembly
      .join(in_reads)
      .set { links_in }
    LINKS(links_in)

    LINKS
      .out
      .scaffolds
      .set { scaffolds }
    MAP_TO_ASSEMBLY(in_reads, scaffolds)

    RUN_QUAST(scaffolds, inputs, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)

    RUN_BUSCO(scaffolds)

    YAK_QC(scaffolds, shortread_kmers)

    if(params.lift_annotations) {
      RUN_LIFTOFF(LINKS.out.scaffolds, inputs)
    }
  emit:
     scaffolds
}

workflow RUN_LONGSTITCH {
  take:
    inputs
    in_reads
    assembly
    references
    ch_aln_to_ref
    shortreads
  
  main:
    assembly
      .join(in_reads)
      .set { longstitch_in }
    LONGSTITCH(longstitch_in)

    LONGSTITCH
      .out
      .ntlLinks_arks_scaffolds
      .set { scaffolds }

    MAP_TO_ASSEMBLY(in_reads, scaffolds)

    RUN_QUAST(scaffolds, inputs, ch_aln_to_ref, MAP_TO_ASSEMBLY.out.aln_to_assembly_bam)

    RUN_BUSCO(scaffolds)

    YAK_QC(scaffolds, shortreads)

    if(params.lift_annotations) {
      RUN_LIFTOFF(LONGSTITCH.out.ntlLinks_arks_scaffolds, inputs)
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
    assembly
    inputs
    aln_to_ref
    aln_to_assembly

  main:
    /* prepare for quast:
     * This makes use of the input channel to obtain the reference and reference annotations
     * See quast module for details
     */
    inputs
      .map { row -> [row.sample, row.ref_fasta, row.ref_gff] }
      .set { inputs_references }

    assembly
      .join(inputs_references)
      .join(aln_to_ref)
      .join(aln_to_assembly)
      .set { quast_in }
    /*
     * Run QUAST
     */
    QUAST(quast_in, use_gff = params.use_ref, use_fasta = false)
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
    BUSCO(assembly, params.busco_lineage, params.busoc_db)
}

/*
 * YAK
 ===========================================
 */


 workflow YAK_QC {
  take:
    assembly
    kmer_shortreads

  main:
    KMER_ASSEMBLY(assembly)
    KMER_ASSEMBLY_HIST()
    if(short_reads) {
      KMER_ASSEMBLY
      .out
      .join(kmer_shortreads)
      .set{ yak_qv_in }
      KMER_ASSEMBLY_QV(yak_qv_in)
    }
 }



/*
 ===========================================
 * ANNOTATIONS
 ===========================================

 * LIFTOFF
 ===========================================
 * Run LIFTOFF to lift annotations from Col-CEN
 */

workflow RUN_LIFTOFF {
  take:
    assembly
    inputs
  
  main:
    assembly
      .join(
        inputs
        .map { row -> [row.sample, row.ref_fasta, row.ref_gff] } )
      .set { liftoff_in }

    LIFTOFF(liftoff_in)
    
    LIFTOFF
      .out
      .set { lifted_annotations }

  emit:
    lifted_annotations
}

 /*
 ====================================================
 ====================================================
                 MAIN WORKFLOW
 ====================================================
 * Collect fastq files
 * PORECHOP
 * NANOQ
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

workflow ASSEMBLE {
 /*
 Define channels
 */

  Channel.empty().set { ch_input }
  Channel.empty().set { ch_refs }
  Channel.empty().set { ch_ref_bam }
  Channel.empty().set { ch_assembly }
  Channel.empty().set { ch_assembly_bam }
  Channel.empty().set { ch_assembly_bam_bai }
  Channel.empty().set { ch_medaka_in }
  Channel.empty().set { ch_polished_genome }
  Channel.empty().set { ch_short_reads }
  Channel.empty().set { sr_kmers }

  /*
  Check samplesheet
  */

  if(params.samplesheet) {
    Channel.fromPath(params.samplesheet) 
           .splitCsv(header:true) 
           .set { ch_input }
    if(params.use_ref) {
      ch_input
        .map { row -> [row.sample, row.ref_fasta] }
        .set { ch_refs }
    }
  } else {
    exit 1, 'Input samplesheet not specified!'
  }

  /*
  Prepare reads
  */
  if(!params.hifi_only) {
    PREPARE_ONT(ch_input)
    JELLYFISH(PREPARE_ONT.out.trimmed, PREPARE_ONT.out.trimmed_med_len)
    PREPARE_ONT
      .out
      .trimmed
      .set { ch_reads }
    ch_ont_reads = ch_reads
    KMER_ONT(ch_ont_reads)
    KMER_ONT
      .out
      .set { ont_kmers }
    ONT_KV(ont_kmers)
    KMER_ONT_HIST(ont_kmers)
  } 
  
  if(short_reads) {
    ch_input
      .map { create_shortread_channel(it) }
      .set { ch_short }
    KMER_SHORTREADS(ch_short)
    KMER_SHORTREADS
      .out
      .set { sr_kmers }
    KMER_SR_HIST(sr_kmers)
  }

  if(params.hifi_only) {
    ch_input
      .map { it -> [it.sample, it.hifireads] }
      .set { ch_reads }
  }
   
  /*
  Assemble with flye, unless --hifi_ont or --hifi_only is set
  */
  if (!params.hifi_ont && !params.hifi_only) {
    ASSEMBLE_ONT(PREPARE_ONT.out.trimmed, ch_input, JELLYFISH.out.hap_len, sr_kmers)
    ASSEMBLE_ONT
      .out
      .ch_assembly
      .set { ch_polished_genome }
    if (params.use_ref) {
      ASSEMBLE_ONT
        .out
        .ch_ref_bam
        .set { ch_ref_bam }
  }
  }
 
  /*
  Optional HiFi assembly
  */

  if(params.hifi) {
    ch_input
      .map { it -> [it.sample, it.hifireads] }
      .set { ch_hifireads }

    if(params.lima) {

      if(is.null(params.pacbio_primers)) error 'Trimming with lima requires a file containing primers (--pacbio_primers)'

      LIMA(ch_hifi_reads, params.pacbio_primers)
      TO_FASTQ(LIMA.out.bam)
      TO_FASTQ
        .out
        .set { ch_hifireads }

    }
    // kmers from (trimmed) reads
    KMER_HIFI(ch_hifireads)
    KMER_HIFI
      .out
      .set { hifi_kmers } 
    HIFI_KV(hifi_kmers)
    KMER_HIFI_HIST(hifi_kmers)

    // If there is a combined assembly for ONT and hifi keep the ont reads.
    if(params.hifi_ont) {
      ch_hifireads
        .join(ch_ont_reads)
        .set { ch_hifireads }
    }

    ASSEMBLE_HIFI(ch_hifireads, ch_input, sr_kmers)
    ASSEMBLE_HIFI
      .out
      .ch_assembly
      .set { ch_hifi_assembly }

    if (params.hifi_ont || params.hifi_only) {
      ch_hifi_assembly 
        .set { ch_polished_genome }
    }
  }
  
  /*
  Polishing with medaka; ONT only
  */

  if(params.polish_medaka) {
    
    if(params.hifi_ont) error 'Medaka should not be used on ONT-HiFi hybrid assemblies'
    if(params.hifi_only) error 'Medaka should not be used on HiFi assemblies'

    POLISH_MEDAKA(ch_input, PREPARE_ONT.out.trimmed, ch_polished_genome, ch_ref_bam, sr_kmers)

    POLISH_MEDAKA
      .out
      .set { ch_polished_genome }
  }

  /*
  Polishing with short reads using pulon
  */

  if(params.polish_pilon) {
    POLISH_PILON(ch_input, ch_reads, ch_polished_genome, ch_ref_bam, sr_kmers)
    POLISH_PILON
      .out
      .pilon_improved
      .set { ch_polished_genome }
  } 

  /*
  Scaffolding
  */

  if(params.scaffold_ragtag) {
    // If there are HiFi and ONT reads assembled individually, ONT will be scaffolded onto HiFi
    if(params.hifi && !params.hifi_ont && !params.hifi_only) {
      RUN_RAGTAG(ch_input, ch_reads, ch_polished_genome, ch_hifi_assembly, ch_ref_bam, sr_kmers)
    } else {
      // In all other casas, the assembly will be scaffolded onto the reference
      RUN_RAGTAG(ch_input, ch_reads, ch_polished_genome, ch_refs, ch_ref_bam, sr_kmers)
    }
  }
  if(params.scaffold_links) {
    RUN_LINKS(ch_input, ch_reads, ch_polished_genome, ch_refs, ch_ref_bam, sr_kmers)
  }
  if(params.scaffold_longstitch) {
    RUN_LONGSTITCH(ch_input, ch_reads, ch_polished_genome, ch_refs, ch_ref_bam, sr_kmers)
  }
  
} 