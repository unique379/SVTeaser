# SVTeaser

SV simulation for rapid benchmarking

[Previous Work](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0803-1)

[Hackathon Schedule](https://docs.google.com/document/d/1ychEMq4vXWtMQRGJD4re5ZzEyIDpb5o_cSBy3CPv2hg/edit#heading=h.5g50ovsn2k70)


## Goals

Make a tool that performs SV and read simulation to create inputs for benchmarking an SV caller. Create an evaluation/reporting procedure of the SV callers’ performance.


## Overview Diagram


![](SVTeaser_Workflow.jpg)


## Working Notes/Documentation

[Here](https://docs.google.com/document/d/1AQxiYEbBhN0-HCAOsrqHZxvsh4ZIFxxeVoJGxApmG-U/edit#)

# File structure diagram 
#### _Define paths, variable names, etc_

## Installation


- Download and install [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR.git)
- Put the `SURVIVOR` executable into your environment's PATH
- Install [vcftools](https://vcftools.github.io/index.html)
- Ensure `vcftools` is in your environment's PATH

## Quick Start

```
usage: svteaser [-h] CMD ...

SVTeaser v0.0.1 - SV simulation for rapid benchmarking

    CMDs:
        sim_sv          Simulate SVs
        surv_sim        Simulate SVs with SURVIVOR
        surv_vcf_fmt    Correct a SURVIVOR simSV vcf
        sim_reads       Run read simulators

positional arguments:
  CMD         Command to execute
  OPTIONS     Options to pass to the command

optional arguments:
  -h, --help  show this help message and exit
```

1. Create a SVTeaser working directory (`output.svt`) by simulating SVs over a reference
2. Simulate reads over the altered reference and place them in the `output.svt` directory
3. Call SVs over the reads (`output.svt/read1.fastq output.svt/read2.fastq`) with your favorite SV caller
4. Run `truvari bench` with the `--base output.svt/simulated.sv.vcf.gz` and `--comp your_calls.vcf.gz`
5. Open the `notebooks/SVTeaser.ipynb` and point to your `output.svt` directory
