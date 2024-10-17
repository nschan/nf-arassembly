include { RUN_PILON } from '../run_pilon/main'
include { MAP_SR } from '../../../mapping/map_sr/main'
include { MAP_TO_ASSEMBLY } from '../../../mapping/map_to_assembly/main'
include { RUN_BUSCO } from '../../../qc/busco/run_busco/main'
include { RUN_QUAST } from '../../../qc/quast/run_quast/main'
include { YAK_QC } from '../../../qc/yak_qc/main'
include { RUN_LIFTOFF } from '../../../liftoff/run_liftoff/main'


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