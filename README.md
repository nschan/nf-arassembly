# nf-arassembly

Assembly pipeline for arabidopsis genomes from nanopore sequencing.

Genome assembly using flye, qc using quast (also works for other species), polish with pilon and / or medaka, scaffold using ragtag, SLR, LINKS or longstitch.

# Parameters

See also [schema.md](schema.md)


| Parameter | Effect |
| --- | --- |
| `--samplesheet` | Path to samplesheet |
| `--collect` | Are the provided reads a folder (true) or a single fq files (default, false) |
| `--flye_mode` | The mode to be used by flye; default: "--nano-hq" |
| `--flye_args` | Arguments to be passed to flye, default: `none`. Example: `--flye_args '--genome-size 130g --asm-coverage 50'` |
| `--medaka_model` | Model used by medaka, default: 'r1041_e82_400bps_hac@v4.2.0:consesus' |
| `--polish_pilon` | Polish with short reads using pilon?; default: false |
| `--busco_db` | Path to local busco db?; default: `/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/software/busco_db` |
| `--skip_flye` | Skip assembly with flye?, requires different samplesheet (!); default: false |
| `--skip_alignments` | Skip alignments? requires different samplesheet (!); default: false |
| `--scaffold_ragtag` | Scaffolding with ragtag?; default: false |
| `--scaffold_links` | Scaffolding with LINKS?; default: false |
| `--scaffold_slr` | Scaffolding with SLR?; default: false |
| `--scaffold_longstitch` | Scaffolding with longstitch?; default: false |
| `--lift_annotations` | Lift annotations from reference?; default: true |
| `--out` | Results directory, default: './results'` |

# Procedure

This pipeline will in a first step extract all fastq.gz files in the readpath folder into a single fastq file. This can be skipped using `--collect false`.

The pipeline then assembles reads into an assembly using flye. 
Then minimap2 aligns reads to the new assembly, and the reference.
These alignments, together with the assembly and reference genome & annotation will then be used as inputs for QUAST.
QUAST will run with the following additional parameters:

```
        --eukaryote \\
        --glimmer \\
        --conserved-genes-finding \\
        --ref-bam ${ref_bam} \\
        --bam ${bam} 
```

Subsequently, the assembly will be polished using first pilon and then medaka, QUAST will again be used to assess the polished genomes.
If --skip_pilon is used the genome will only be polished using medaka


# Usage

## Full Pipeline

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
`--polish_pilon true`

This requires additional information in the samplesheet, shortread_F, shortread_R and paired

```
sample,readpath,assembly,ref_fasta,ref_gff,assembly_bam,assembly_bai,ref_bam,shortread_F,shortread_R,paired
sampleName,reads,assembly.fasta.gz,reference.fasta,reference.gff,reads_on_assembly.bam,reads_on_assembly.bai,reads_on_reference.bam,short_F1.fastq,short_F2.fastq,true
```

In a case where only single-reads are available, shortread_R should be empty, and paired should be false

## Scaffolding

`ragtag`, `LINKS`, `SLR` and `longstitch` can be used for scaffolding.
The `SLR` container segfaults with illegal instructions, I assume that the container is somehow built wrong. 

## Using liftoff

If `lift_annotations` is used, the annotations from the reference genome will be mapped to assemblies and scaffolds using liftoff. This will happen at each step of the pipeline where a new genome version is created, i.e. after assembly, after polishing and after scaffolding.

# Problems

Presumably, something like this will happen:

```
Error executing process > 'ALIGN_TO_REF (TestSample)'

Caused by:
  Charliecloud failed to pull image
  command: ch-image pull -s /dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/charliecloud/storage/0.30 gitlab.lrz.de:5005/beckerlab/container-playground/minimap2-samtools:74bccfe8 > /dev/null
  status : 1
  message:
    pulling image:    gitlab.lrz.de:5005/beckerlab/container-playground/minimap2-samtools:74bccfe8
    requesting arch:  amd64
    error: unauthorized or not in registry: gitlab.lrz.de:5005/beckerlab/container-playground/minimap2-samtools:74bccfe8
    hint: if your registry needs authentication, use --auth
```

Follow the "hint" and just copy the `ch-image pull -s ...` line (excluding `> /dev/null`) and append `--auth` to pull the container, authenticate with lrz credentials, then retry.
