# nf-arassembly

Arabidopsis assembly using flye, qc using quast (also works for other species).

# Procedure

This pipeline will in a first step (COLLECT_READS) extract all fastq.gz files in the readpath folder into a single fastq file.
The pipeline then assembles using flye and uses minimap2 to align reads to the new assembly, and the reference.
These alignments, together with the assembly and reference genome & annotation will then be used as inputs for QUAST.
QUAST will run with the following additional parameters:

```
        --eukaryote \\
        --gene-finding \\
        --conserved-genes-finding \\
        --ref-bam ${ref_bam} \\
        --bam ${bam} 
```

# Usage

## With a samplesheet

I assume that creating a samplesheet is the easiest way to do this.

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

## Parameters

`--out`: output folder, typically './results'. Needs to be provided.

`--samplesheet`: path to the samplesheet (required)

`--collect`: should reads be collected into a single fastq? Default: true

`--flye_mode`: changes mode for flye, default is "--nano-hq".
Valid options are: "--pacbio-raw", "--pacbio-corr", "--pacbio-hifi", "--nano-raw", "--nano-corr", "--nano-hq"

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
