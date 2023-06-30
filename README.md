# nf-arassembly

Genome assembly using flye, qc using quast (also works for other species), polish with pilon and / or medaka.

# Parameters

`--samplesheet`            Path to samplesheet

`--collect`                    Are the provided reads a folder (default, true) or a single fq files (false)

`--skip_flye`                Should the flye assembly step be skipped (default false), requires different samplesheet (!)

`--flye_mode`                The mode to be used by flye; default: "--nano-hq"

`--medaka_model`          Model used by medaka, default: 'r1041_e82_400bps_hac_v4.2.0'

`--skip_pilon`              Should pilon be skipped; default: false 

`--skip_alignments`   Should the alignments be skipped, requires different samplesheet (!), default: false

`--out`                             Results directory, default: './results'`

# Procedure

This pipeline will in a first step (COLLECT_READS) extract all fastq.gz files in the readpath folder into a single fastq file.
The pipeline then assembles using flye and uses minimap2 to align reads to the new assembly, and the reference.
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
