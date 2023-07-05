include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LIFTOFF {
  tag "$meta"
  label 'process_high'
  
  conda "bioconda::liftoff=1.6.4"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/liftoff:1.6.3--pyhdfd78af_0' :
          'quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(assembly), path(reference_fasta), path(reference_gff)
  
  output:
      tuple val(meta), path("${assembly}_liftoff.gff"), emit: lifted_annotations

  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  liftoff \\
    -g ${reference_gff} \\
    ${assembly} \\
    ${reference_fasta} \\
    -o ${assembly}_liftoff.gff
  """
}