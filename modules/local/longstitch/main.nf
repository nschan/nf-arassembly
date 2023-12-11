include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LONGSTITCH {
  tag "$meta"
  label 'process_high'
  
  publishDir "${params.out}",
      mode: params.publish_dir_mode,
      saveAs: { filename -> saveFiles(filename:filename,
                                      options:params.options, 
                                      publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                      publish_id:meta) }
  input:
      tuple val(meta), path(assembly), path(reads)

  output:
      tuple val(meta), path("soft-links/*longstitch-scaffolds.fa"), emit: scaffolds
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  longstitch tigmint-ntLink-arks draft=${assembly} reads=${reads} t=${task.cpus} G=135e6 out_prefix=${meta}
  """
}