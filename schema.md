# . pipeline parameters



## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |

## Other parameters

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `samplesheet` | samplesheet | `string` |  | True |  |
| `enable_conda` | use conda (not supported) | `boolean` |  |  |  |
| `collect` | combine reads from multiple fastq files into one? | `boolean` |  |  |  |
| `skip_flye` | skip flye assembly | `boolean` |  |  |  |
| `skip_alignments` | skip alignment to initial assembly and reference | `boolean` |  |  |  |
| `flye_mode` | flye mode to use | `string` | '--nano-hq' |  |  |
| `flye_args` | extra args for flye | `string` | '' |  |  |
| `polish_pilon` | use pilon to polish with short reads | `boolean` |  |  |  |
| `medaka_model` | medaka model to use | `string` | 'r1041_e82_400bps_hac@v4.2.0:consesus' |  |  |
| `lift_annotations` | run liftoff to lift over annotations | `boolean` |  |  |  |
| `out` | results directory | `string` | './out' |  |  |
| `scaffold_ragtag` | Scaffold with ragtag | `boolean` |  |  |  |
| `scaffold_links` | Scaffold with links | `boolean` |  |  |  |
| `scaffold_longstitch` | Scaffold with longstitch | `boolean` |  |  |  |
| `busoc_db` | path to busco db | `string` | /dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/software/busco_db |  |  |
| `porechop` | run porechop | `boolean` |  |  |  |
| `use_ref` | Use a reference genome | `boolean` | True |  |  |
| `busco_lineage` | busco lineage to use | `string` | "brassicales_odb10" |  |  |
