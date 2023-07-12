include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process NANOQ {
    tag "$meta"
    label 'process_low'
       publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta) }

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*_report.json"), emit: report
        tuple val(meta), path("*_stats.json"), emit: stats

    script:
        """
        nanoq -i ${reads} -j -r ${meta}_report.json -s -H -vvv > ${meta}_stats.json
        """
}