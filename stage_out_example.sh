#!/bin/bash

# This file is used to stage out the relevant output from nf-arassembly into a storage directory
# Main purpose is to skip aligment files and other big files.
# It also makes the names a bit easier to handle.
# Niklas Schandry; 2024

accessions=(
    Bor-4
    Borky-1
    Bozen-1_2
    Faeb-2
    IP_Hue-3
    IP_Var-0
    Istisu-1
    Shigu-2
    T960
    Tac-0
    Ts-1
    Ukidi-07
    Kidr
    Vimmerby
    T840
    Castelfed-1-199
    Hadd-3
    Col-0
    Paw-26
    Rak-2
    Kastel-1
    Tu-B1-2
    WI-0
    DraIV3-7
    Giffo-1
    Dra2-1
    Borky-1-new
)

stages=(
    assembly
    polish_medaka
    run_links
    run_longstitch
    run_ragtag
)

dest_base="/dss/dsslegfs01/pn73so/pn73so-dss-0000/becker_common/Duncan/assemblies/"

results_base="/dss/lxclscratch/0A/ra34bin/assembly_chop/results/arassembly/"
for acc in "${accessions[@]}"; do
    for stage in "${stages[@]}"; do

        if [[ "$stage" == "assembly" ]]; then
            fasta_dir="${results_base}/${stage}/flye"
        elif [[ "$stage" == "polish_medaka" ]]; then
            fasta_dir="${results_base}/${stage}/run_medaka/medaka"
        elif [[ "$stage" == "run_links" ]]; then
            fasta_dir="${results_base}/${stage}/links"
        elif [[ "$stage" == "run_longstitch" ]]; then
            fasta_dir="${results_base}/${stage}/longstitch"
        elif [[ "$stage" == "run_ragtag" ]]; then
            fasta_dir="${results_base}/${stage}/ragtag_scaffold/${acc}*/"
        fi

        echo "Moving $acc stage $stage to $dest_base/$stage"
        echo "  Moving fasta files"
        mkdir -p ${dest_base}/${stage}/fasta
        cp -n ${fasta_dir}/${acc}*.fa* ${dest_base}/${stage}/fasta/ 2>/dev/null || :

        echo "  Moving annotations"
        mkdir -p ${dest_base}/${stage}/annotations/liftoff
        cp -n ${results_base}/${stage}/run_liftoff/liftoff/${acc}*.gff ${dest_base}/${stage}/annotations/liftoff/ 2>/dev/null || :

        echo "  Moving busco"
        mkdir -p ${dest_base}/${stage}/busco
        cp -n ${results_base}/${stage}/run_busco/busco/short_summary*${acc}*.fa*.txt ${dest_base}/${stage}/busco/ 2>/dev/null || :

        echo "  Moving quast"
        mkdir -p ${dest_base}/${stage}/quast
        cp -n ${results_base}/${stage}/run_quast/quast/${acc}/report.tsv ${dest_base}/${stage}/quast/${acc}_report.tsv 2>/dev/null || :
        cp -n ${results_base}/${stage}/run_quast/quast/${acc}/report.pdf ${dest_base}/${stage}/quast/${acc}_report.pdf 2>/dev/null || :
    done
done
