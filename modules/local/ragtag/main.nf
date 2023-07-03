include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RAGTAG_SCAFFOLD {
  tag "$meta"
  label 'process_high'
  
  conda "bioconda::ragtag=2.1.0"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/ragtag:2.1.0--pyhb7b1952_0' :
          'quay.io/biocontainers/ragtag:2.1.0--pyhb7b1952_0' }"

  publishDir "${params.out}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

  input:
      tuple val(meta), path(assembly), path(reference)
  
  output:
      tuple val(meta), path("${assembly}_corrected_${reference}/*.fasta"), emit: corrected_assembly
      tuple val(meta), path("${assembly}_corrected_${reference}/*.agp"),   emit: corrected_agp
      tuple val(meta), path("${assembly}_corrected_${reference}/*.stats"), emit: corrected_stats
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  ragtag.py scaffold ${assembly} ${reference} \\
    -o "${assembly}_corrected_${reference}" \\
    -t $task.cpus \\
    -f 15000 \\
    --remove-small

  mv ${assembly}_corrected_${reference}/ragtag.scaffold.fasta ${assembly}_corrected_${reference}/${assembly}_corrected_${reference}.fasta
  """
}