process MERYL_COUNT {
    tag "$meta"
    label 'process_high'
    publishDir(
      path: { "${params.out}/merqury/meryl/" }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
    tuple val(meta), path(reads)
    val kvalue

    output:
    tuple val(meta), path("*.meryl")    , emit: meryl_db

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    for READ in $reads; do
        meryl count \\
            k=$kvalue \\
            threads=$task.cpus \\
            memory=${task.memory.toGiga()} \\
            $args \\
            $reads \\
            output read.\${READ%.f*}.meryl
    done
    """
}