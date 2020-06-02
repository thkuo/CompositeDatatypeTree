#!/usr/bin/env bash

source activate snakemake_env
snakemake --configfile ./config.yml --snakefile ./gpa.smk ../../results.v5/gpa/gpa.var.aln
