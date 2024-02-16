
config="/shared/home/jpeter/Scripts/FLEPseq2/Genes_list/yeast_flepseq/config_yeast.yaml"
#config="../config/config_FLEPseq_SpikeIn.yaml"

#snakemake --version
snakemake --profile slurm --use-conda --configfile $config -s SnakefileFLEPseq2
snakemake --profile slurm --use-conda --configfile $config -s SnakefileDownstream

