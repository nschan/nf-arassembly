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
      tuple val(meta), path("${assembly}.tigmint-ntLink-arks.longstitch-scaffolds.fa"), emit: ntlLinks_arks_scaffolds
      tuple val(meta), path("${assembly}.tigmint-ntLink.longstitch-scaffolds.fa"), emit: ntlLinks_scaffolds
  
  script:
      def prefix = task.ext.prefix ?: "${meta}"
  """
  if [[ ${assembly} == *.gz ]]; then
    zcat ${assembly} > assembly.fasta
  fi

  if [[ ${assembly} == *.fa || ${assembly} == *.fasta ]]; then
    cp ${assembly} assembly.fasta
  fi
  ln -s assembly.fasta assembly.fa

  if [[ ${reads} == *.fa || ${reads} == *.fasta ]]; then
    gzip -c ${reads} > ${reads}.gz
    ln -s ${reads}.gz reads.fq.gz
  fi

  if [[ ${reads} == *.gz ]]; then
    cp ${reads} reads.fq.gz
  fi

  longstitch tigmint-ntLink-arks draft=assembly reads=reads t=${task.cpus} G=135e6 out_prefix=${assembly}
  cp assembly.k32.w100.tigmint-ntLink-arks.longstitch-scaffolds.fa ${assembly}.tigmint-ntLink-arks.longstitch-scaffolds.fa
  cp assembly.k32.w100.tigmint-ntLink.longstitch-scaffolds.fa ${assembly}.tigmint-ntLink.longstitch-scaffolds.fa
  """
}