include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LINKS {
  tag "$meta"
  label 'process_high'
  
  container "quay.io/biocontainers/links:2.0.1--h9f5acd7_3"

  publishDir "${params.out}",
      mode: params.publish_dir_mode,
      saveAs: { filename -> saveFiles(filename:filename,
                                      options:params.options, 
                                      publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                      publish_id:meta) }
  input:
      tuple val(meta), path(assembly), path(reads)

  output:
      tuple val(meta), path("*.scaffolds.fa"), emit: scaffolds
      tuple val(meta), path("*.scaffolds"), emit: scaffold_csv
      tuple val(meta), path("*.gv"), emit: graph
      tuple val(meta), path("*.log"), emit: log
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  LINKS -f ${assembly} -s ${reads} -j 4 -b ${meta}_links
  """
}