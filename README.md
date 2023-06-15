# nf-remap

Map your long reads to assemblies.

# Why

You made assemblies from long reads, now you want to map the reads to the assemblies.

# How


## With a samplesheet

I assume that creating a samplesheet is the easiest way to do this.

The samplesheet _must_ adhere to this format, including the header row. Please note the absence of spaces after the commas:

```
sample,reads,ref
Name,path/to/Sample_Reads.fastq,path/to/Reference.fasta
```

To run the pipeline with a samplesheet on biohpc_gen:
```
nextflow run reseq-nf/  --samplesheet 'path/to/sample_sheet.csv' \
                        --out './results' \
                        -profile charliecloud,biohpc_gen
```

## Working example

Samplesheet looks like this and is saved as samplesheet.csv:

```
sample,reads,ref
barcode49,/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/read_data/A_thaliana_ONT/single_fq/barcode49.fastq,/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/Duncan/assemblies/flye/barcode49/OUT/assembly.fasta

```
Nextflow call like this:

```
nextflow run ~/nf-remap --samplesheet samplesheet.csv -profile charliecloud,biohpc_gen
```

Presumably, something like this will happen:

```
Error executing process > 'ALIGN (barcode49)'

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

Follow the "hint" and just copy the `ch-image pull -s ...` line and append `--auth` to pull the container, authenticate with lrz credentials, then retry.
