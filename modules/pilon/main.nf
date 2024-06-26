process PILON {
    tag "$meta"
    label 'process_medium'
    conda "bioconda::pilon=1.24"
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 

    input:
    tuple val(meta), path(fasta), path(bam), path(bai)
    val pilon_mode

    output:
    tuple val(meta), path("*.fasta.gz") , emit: improved_assembly
    tuple val(meta), path("*.vcf")      , emit: vcf               , optional : true
    tuple val(meta), path("*.change")   , emit: change_record     , optional : true
    tuple val(meta), path("*.bed")      , emit: tracks_bed        , optional : true
    tuple val(meta), path("*.wig")      , emit: tracks_wig        , optional : true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def valid_mode = ["frags", "jumps", "unpaired", "bam"]
    if ( !valid_mode.contains(pilon_mode) )  { error "Unrecognised mode to run Pilon. Options: ${valid_mode.join(', ')}" }
    """
    pilon \\
        --genome $fasta \\
        --output ${meta}_pilon \\
        --threads $task.cpus \\
        $args \\
        --$pilon_mode $bam \\
        -Xmx90G

    gzip -n ${prefix}_pilon.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pilon: \$(echo \$(pilon --version) | sed 's/^.*version //; s/ .*\$//' )
    """
}