# nf-arassembly

Assembly pipeline for arabidopsis genomes from nanopore sequencing. Should also work for other species.

Genome assembly using porechop, flye, qc using quast, polish with pilon and / or medaka, scaffold using ragtag, LINKS or longstitch.

# Parameters

See also [schema.md](schema.md)


| Parameter | Effect |
| --- | --- |
| `--samplesheet` | Path to samplesheet |
| `--collect` | Are the provided reads a folder (`true`) or a single fq files (default: `false` ) |
| `--porechop` | Run porechop? (default: `false`) |
| `--use_ref` | Use a refence genome? (default: `true`) |
| `--flye_mode` | The mode to be used by flye; default: `"--nano-hq"` |
| `--flye_args` | Arguments to be passed to flye, default: `none`. Example: `--flye_args '--genome-size 130g --asm-coverage 50'` |
| `--medaka_model` | Model used by medaka, default: `'r1041_e82_400bps_hac@v4.2.0:consesus'` |
| `--polish_pilon` | Polish with short reads using pilon? Sefault: `false` |
| `--busco_db` | Path to local busco db?; default: `/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/software/busco_db` |
| `--busco_lineage` | Busco lineage to use; default: `brassicales_odb10` |
| `--skip_flye` | Skip assembly with flye?, requires different samplesheet (!); Default: `false` |
| `--skip_alignments` | Skip alignments? requires different samplesheet (!); Default: `false` |
| `--scaffold_ragtag` | Scaffolding with ragtag? Default: `false` |
| `--scaffold_links` | Scaffolding with LINKS? Default: `false` |
| `--scaffold_slr` | Scaffolding with SLR? Default: `false` |
| `--scaffold_longstitch` | Scaffolding with longstitch? Default: `false` |
| `--lift_annotations` | Lift annotations from reference using liftoff? Default: `true` |
| `--out` | Results directory, default: `'./results'` |

# Procedure

  * Extract all fastq.gz files in the readpath folder into a single fastq file. By default this is skipped, enable with `--collect`.
  * Barcodes and adaptors will be removed using porechop. By default this is skipped, enable with `--porechop`.
  * Read QC is done via nanoq.
  * Assemblies are performed with flye.
  * Polishing is done using medaka, and scaffolding via links, longstitch and / or ragtag. 
  * Optional short-read polishing can be done using pilon. By default this is not done, enable with `--polish_pilon`, requires different samplesheet with shortreads.
  * Annotations are created using liftoff. 
  * Quality of each stage is assessed using QUAST and BUSCO (standalone).

# Graph

![Tubemap](nf-arassembly.tubes.png)

# Additional information

QUAST will run with the following additional parameters:

```
        --eukaryote \\
        --glimmer \\
        --conserved-genes-finding \\
```

# Usage

## Standard Pipeline

The samplesheet _must_ adhere to this format, including the header row. Please note the absence of spaces after the commas:

```
sample,readpath,ref_fasta,ref_gff
sampleName,path/to/reads,path/to/reference.fasta,path/to/reference.gff
```

To run the pipeline with a samplesheet on biohpc_gen:
```
nextflow run nf-arassembly --samplesheet 'path/to/sample_sheet.csv' \
                           --out './results' \
                           -profile charliecloud,biohpc_gen
```

## No refence genome

If there is no reference genome available use `--use_ref false` to disable the reference genome.
Liftoff should not be used without a reference, QUAST will no longer compare to reference. 

## Skipping Flye

In case you already have an assembly and would only like to check it with QUAST and polish use
`--skip_flye true`

This mode requires a different samplesheet:

```
sample,readpath,assembly,ref_fasta,ref_gff
sampleName,path/to/reads,assembly.fasta.gz,reference.fasta,reference.gff
```

When skipping flye the original reads will be mapped to the assembly and the reference genome.

## Skipping Flye and mappings

In case you have an assembly and have already mapped your reads to the assembly and the reference genome you can use
`--skip_flye true --skip_alignments true`

This mode requires a different samplesheet:

```
sample,readpath,assembly,ref_fasta,ref_gff,assembly_bam,assembly_bai,ref_bam
sampleName,reads,assembly.fasta.gz,reference.fasta,reference.gff,reads_on_assembly.bam,reads_on_assembly.bai,reads_on_reference.bam
```

## Polishing with pilon

The assemblies can optionally be polished using available short-reads using pilon.
`--polish_pilon`

This requires additional information in the samplesheet: `shortread_F`, `shortread_R` and `paired`:

```
sample,readpath,ref_fasta,ref_gff,shortread_F,shortread_R,paired
sampleName,reads,assembly.fasta.gz,reference.fasta,reference.gff,short_F1.fastq,short_F2.fastq,true
```

In a case where only single-reads are available, `shortread_R` should be empty, and `paired` should be false

## Scaffolding

`ragtag`, `LINKS`, and `longstitch` can be used for scaffolding.

## Using liftoff

If `lift_annotations` is used, the annotations from the reference genome will be mapped to assemblies and scaffolds using liftoff.
This will happen at each step of the pipeline where a new genome fasta is created, i.e. after assembly, after polishing and after scaffolding.