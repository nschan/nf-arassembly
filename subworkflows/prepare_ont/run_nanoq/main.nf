include { NANOQ } from '../../../modules/nanoq/main'

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
    .median_length
    .set { median_length }

  emit:
    report
    stats
    median_length
 }