
config="/shared/home/jpeter/Scripts/FLEPseq2/Genes_list/human_flepseq/config_human.yaml"
#config="../config/config_FLEPseq_SpikeIn.yaml"

#snakemake --version
snakemake --profile slurm --use-conda --configfile $config -s SnakefileFLEPseq2
snakemake --profile slurm --use-conda --configfile $config -s SnakefileDownstream

