# FLEPseq2
## TO DO: Complete this file. 

## Description

Snakemake pipelines for polyA detection via FLEPseq2

These snakemakes regroup the major steps described in the FLEPSeq2 github repository (https://github.com/ZhaiLab-SUSTech/FLEPSeq). 
We added an extra step (extract_tails.py) to extract the additionnal tail after the polyA, and analyse its composition and a downstream analysis for graph generation.

Some minor changes have been done to the original FLEPSeq2 code:
- Handling directly FASTQ files without needing to convert them to FASTA

Steps of the workflow
```mermaid
graph TD
    A(.nanopore.fast5) --> |Basecalling with Guppy| B(.fastq)
    B --> |Mapping with minimap2| C(.bam)
    C --> |adapterFinder.py| D(.adapter_results.txt)
    D --> |PolyA Caller| E(results)
    D --> |Downstream Analysis|E(results)
    D --> |3' terminal non-adenosine analysis|E(results)
    E --> |plot_tail.R| F(result.merged.parts.csv)

```
## Getting started

Clone the repository

```
git clone https://github.com/jackson-peter/FLEPseq2.git
```

A configuration file (.yaml) is required to run the workflow. An example is included in the FLEPseq2/config folder. You can modify it to suit your data.
```
basecalled_dir: "/ssd_workspace/jpeter/ssData/Guppy_basecalling/RUN11_F/workspace/"
outdir: "/home/jpeter/DATA/FLEPseq/RUN11_F/"
barcode_corr: "/home/jpeter/DATA/FLEPseq/RUN11_F/barcode_correspondance.tsv"

runs_dir: "1_Runs"
mapping_dir: "2_Mapping"
polyA_dir: "3_PolyA"
tail_dir: "4_Tail"

introns_exons: "/home/jpeter/DATA/ReferenceGenomes/Athaliana/TAIR10/exon_intron_pos.repr.bed"
select_introns: "/home/jpeter/DATA/ReferenceGenomes/Athaliana/TAIR10/select_introns.txt"

reference_genome: "/home/jpeter/DATA/ReferenceGenomes/Athaliana/TAIR10/TAIR10_chr_all.fas" # Reference genome in fasta

## MINIMAP MAPPING ADDITIONAL PARAMETERS
minimap2_add_opts: "--secondary=no -G 5000"

```

## Usage

To run the whole pipeline, just execute the runFLEPseq.sh file containing the snakemakes commands
```
conda activate snakemake # or mamba activate snakemake if it is installed
# Go in the 'workflow' folder where the Snakefile is 
cd FLEPseq2/workflow
bash runFLEPseq.sh
```

## Requirements

The only requirement to run the snakefiles is:
- snakemake (>=7)

The other requirements will be handled directly by snakemake with the --use-conda option. This will automatically set up an environement containing all the requirements for the pipeline's execution.



## Usefull links

- Research article describing FLEPSeq2: https://www.nature.com/articles/s41596-021-00581-7
