process COLLECT_READS {
    tag "$meta"
    label 'process_low'

    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    conda "conda-forge::python=3.8.3"
  
  input:
    tuple val(meta), path(read_directory)
  
  output:
    tuple val(meta), path("*.fastq"), emit: combined_reads
  
  script:
    def prefix = task.ext.prefix ?: "${meta}"
      
  """
  gunzip -c ${read_directory}/*.gz > ${prefix}_all_reads.fastq
  """
}