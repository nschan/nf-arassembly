#!/usr/bin/env nextflow

/*
Parameter setup
*/
nextflow.enable.dsl = 2 
params.publish_dir_mode = 'copy'
params.samplesheet = false
params.enable_conda = false
params.collect = false
params.porechop = false
//Jellyfish params
params.jelly_is_reads = true
params.kmer_length = 21
params.read_length = null
params.dump = false

//
params.use_ref = true
params.skip_assembly = false
params.genome_size = null
params.flye_mode = '--nano-hq'
params.flye_args = ''
params.hifi = false
params.lima = false
params.pacbio_primers = null
params.assembler = "flye"
params.qc_reads = null
params.hifiasm_ont = false
params.hifi_only = false
params.hifi_args = ''
params.short_reads = false
params.polish_pilon = false
params.polish_medaka = true
params.medaka_model = 'r1041_e82_400bps_hac_v4.2.0'
params.skip_alignments = false
params.scaffold_ragtag = false
params.scaffold_links = false
params.scaffold_longstitch = false
params.lift_annotations = true
params.busoc_db = ""
params.busco_lineage = "brassicales_odb10"
params.out = './results'

/*
 Print very cool text and parameter info to log. 
*/

log.info """\
======================================================================================================================================================
======================================================================================================================================================
███▄▄▄▄      ▄████████    ▄████████    ▄████████    ▄████████    ▄████████    ▄████████    ▄████████   ▄▄▄▄███▄▄▄▄   ▀█████████▄   ▄█       ▄██   ▄   
███▀▀▀██▄   ███    ███   ███    ███   ███    ███   ███    ███   ███    ███   ███    ███   ███    ███ ▄██▀▀▀███▀▀▀██▄   ███    ███ ███       ███   ██▄ 
███   ███   ███    █▀    ███    ███   ███    ███   ███    ███   ███    █▀    ███    █▀    ███    █▀  ███   ███   ███   ███    ███ ███       ███▄▄▄███ 
███   ███  ▄███▄▄▄       ███    ███  ▄███▄▄▄▄██▀   ███    ███   ███          ███         ▄███▄▄▄     ███   ███   ███  ▄███▄▄▄██▀  ███       ▀▀▀▀▀▀███ 
███   ███ ▀▀███▀▀▀     ▀███████████ ▀▀███▀▀▀▀▀   ▀███████████ ▀███████████ ▀███████████ ▀▀███▀▀▀     ███   ███   ███ ▀▀███▀▀▀██▄  ███       ▄██   ███ 
███   ███   ███          ███    ███ ▀███████████   ███    ███          ███          ███   ███    █▄  ███   ███   ███   ███    ██▄ ███       ███   ███ 
███   ███   ███          ███    ███   ███    ███   ███    ███    ▄█    ███    ▄█    ███   ███    ███ ███   ███   ███   ███    ███ ███▌    ▄ ███   ███ 
 ▀█   █▀    ███          ███    █▀    ███    ███   ███    █▀   ▄████████▀   ▄████████▀    ██████████  ▀█   ███   █▀  ▄█████████▀  █████▄▄██  ▀█████▀  
                                      ███    ███                                                                                  ▀                   
------------------------------------------------------------------------------------------------------------------------------------------------------
Niklas Schandry                                          niklas@bio.lmu.de                                                                    
------------------------------------------------------------------------------------------------------------------------------------------------------
  Results directory  : ${params.out}

  General parameters
     samplesheet     : ${params.samplesheet}
     use reference   : ${params.use_ref}

  ONT preprocessing
     collect         : ${params.collect}
     porechop        : ${params.porechop}

  pacbio preprocessing
    lima             : ${params.lima}
    pacbio primers   : ${params.pacbio_primers}
   
  Assembler          : ${params.assembler}

  Flye assembly
     read_length     : ${params.read_length}
     genome_size     : ${params.genome_size}
     flye_mode       : ${params.flye_mode}

  Hifiasm assembly   : ${params.hifi} 
     Mix HiFi and ONT: ${params.hifiasm_ont}
     Only HiFi       : ${params.hifi_only}
     hifiasm args    : ${params.hifi_args}

  ONT Polishing
     Run Medaka      : ${params.polish_medaka}
     Medaka model    : ${params.medaka_model}

  Short-reads        : ${params.short_reads} 
     Run pilon       : ${params.polish_pilon}

  BUSCO parameters
     busco db        : ${params.busoc_db}
     busco lineage   : ${params.busco_lineage}
     use reference   : ${params.use_ref}

  Scaffolding Tools
     ragtag          : ${params.scaffold_ragtag}
     LINKS           : ${params.scaffold_links}
     longstitch      : ${params.scaffold_longstitch}

  Annotation lift    : ${params.lift_annotations}

  Steps skipped
     flye            : ${params.skip_assembly}
     alignments      : ${params.skip_alignments}
======================================================================================================================================================
======================================================================================================================================================
"""
    .stripIndent(false)

include { ASSEMBLE } from './subworkflows/main'

workflow {
  ASSEMBLE()
}