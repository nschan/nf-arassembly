include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GENOMESCOPE {
    tag "$meta"
    label 'process_medium'
    publishDir "${params.out}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename,
                                        options:params.options, 
                                        publish_dir:"${task.process}".replace(':','/').toLowerCase(), 
                                        publish_id:meta) }
    input:
        tuple val(meta), path(histo), val(kmer_length), val(read_length)

    output:
        tuple val(meta), path("*_genomescope.txt"), emit: summary
        tuple val(meta), path("*_plot.log.png"), emit: plot_log
        tuple val(meta), path("*_plot.png"), emit: plot

    def complement = params.is_reads ? '-C' : ''

    script:
        """
        genomescope.R $histo $kmer_length $read_length genomescope
        mv genomescope/summary.txt ${meta}_genomescope.txt
        mv genomescope/plot.log.png ${meta}_plot.log.png
        mv genomescope/plot.png ${meta}_plot.png
        """
}
