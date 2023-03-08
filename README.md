# FLEPseq2
## TO DO: Complete this file. 

## Description

Snakemake pipeline for polyA detection via FLEPseq2

This snakemake regroups the major steps described in the FLEPSeq2 github repository (https://github.com/ZhaiLab-SUSTech/FLEPSeq). 
We added an extra step (extract_tails.py) to extract the additionnal tail after the polyA, and analyse its composition. 

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

## Usefull links

- Research article describing FLEPSeq2: https://www.nature.com/articles/s41596-021-00581-7
