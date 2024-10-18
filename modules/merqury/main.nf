process MERQURY {
    tag "$meta"
    label 'process_low'
    publishDir(
      path: { "${params.out}/merqury/merqury/" }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
    tuple val(meta), path(meryl_db), path(assembly)

    output:
    tuple val(meta), path("*_only.bed")          , emit: assembly_only_kmers_bed
    tuple val(meta), path("*_only.wig")          , emit: assembly_only_kmers_wig
    tuple val(meta), path("*.completeness.stats"), emit: stats
    tuple val(meta), path("*.dist_only.hist")    , emit: dist_hist
    tuple val(meta), path("*.spectra-cn.fl.png") , emit: spectra_cn_fl_png
    tuple val(meta), path("*.spectra-cn.hist")   , emit: spectra_cn_hist
    tuple val(meta), path("*.spectra-cn.ln.png") , emit: spectra_cn_ln_png
    tuple val(meta), path("*.spectra-cn.st.png") , emit: spectra_cn_st_png
    tuple val(meta), path("*.spectra-asm.fl.png"), emit: spectra_asm_fl_png
    tuple val(meta), path("*.spectra-asm.hist")  , emit: spectra_asm_hist
    tuple val(meta), path("*.spectra-asm.ln.png"), emit: spectra_asm_ln_png
    tuple val(meta), path("*.spectra-asm.st.png"), emit: spectra_asm_st_png
    tuple val(meta), path("${prefix}.qv")        , emit: assembly_qv
    tuple val(meta), path("${prefix}.*.qv")      , emit: scaffold_qv
    tuple val(meta), path("*.hist.ploidy")       , emit: read_ploidy
    tuple val(meta), path("*.hapmers.blob.png")  , emit: hapmers_blob_png           , optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta}"
    def VERSION = 1.3 // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi
    # limit meryl to use the assigned number of cores.
    export OMP_NUM_THREADS=$task.cpus

    merqury.sh \\
        $meryl_db \\
        $assembly \\
        $prefix
    """
}