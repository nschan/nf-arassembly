# . pipeline parameters



## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |

## Other parameters

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `samplesheet` |  | `string` | false |  |  |
| `enable_conda` |  | `string` | false |  |  |
| `collect` |  | `string` | false |  |  |
| `skip_flye` |  | `string` | false |  |  |
| `skip_alignments` |  | `string` | false |  |  |
| `flye_mode` |  | `string` | '--nano-hq' |  |  |
| `flye_args` |  | `string` | '' |  |  |
| `polish_pilon` |  | `string` | false |  |  |
| `medaka_model` |  | `string` | 'r1041_e82_400bps_hac_v4.2.0' |  |  |
| `lift_annotations` |  | `string` | false |  |  |
| `out` |  | `string` | './out' |  |  |
| `scaffold_ragtag` |  | `string` | false |  |  |
| `scaffold_links` |  | `string` | false |  |  |
| `scaffold_slr` |  | `string` | false |  |  |
| `scaffold_longstitch` |  | `string` | false |  |  |
