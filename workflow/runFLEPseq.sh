
config="../config/config_FLEPseq2.yaml"
#config="../config/config_FLEPseq_SpikeIn.yaml"

snakemake --profile slurm --use-conda --configfile $config -s SnakefileFLEPseq2
snakemake --profile slurm --use-conda --configfile $config -s SnakefileDownstream

