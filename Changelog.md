# Changelog

All notable changes to this project will be documented in this file.

The version of the pipeline is indicated in the version.log output file of the pipeline.
You can also get it by executing git describe command with the tags, long & dirty flags.

The pipeline tries to respect semver specs (https://semver.org/). 

```bash
# v0.1.0-1-gdeadbee-dirty
# ^      ^ ^^       ^
# |      | ||       |
# |      | ||       '-- flag indicating if local copy is dirty or not
# |      | |'-- SHA of HEAD (first seven chars)
# |      | '-- "g" is for git
# |      '---- number of commits since last tag
# |
# '--------- last tag
```
## [v1.0.0] 18/04/2023

First stable release of the pipeline. As most of the impacting modifications affect the output of extract_tails.py script, here is the header of this output:

* read_core_id
* chr
* read_exon_total_num
* mRNA
* mRNA_intron_num
* mRNA_start
* mRNA_end
* retention_introns
* primer_type
* polya_start_raw
* polya_end_raw
* polya_start_base
* polya_end_base
* polya_length
* init_polya_start_base
* init_polya_end_base
* init_polya_length
* type
* readname
* sense
* polytail
* additional_tail
* adapter
* dist_adapter
* coords_in_read
* comment
* A_tail_A
* A_tail_T
* A_tail_G
* A_tail_C
* A_tail_pct_A
* A_tail_pct_T
* A_tail_pct_G
* A_tail_pct_C
* add_tail_A
* add_tail_T
* add_tail_G
* add_tail_C
* add_tail_pct_A
* add_tail_pct_T
* add_tail_pct_G
* add_tail_pct_C 

## [v0.2.0] 15/04/2023

### Added
 - a log for the merge_read_info.R to track the different steps (and the number of lines of each step)
 - additional tail detection even when there is no polyA tail (some transcripts like the elongating_3_mapping_low_accuracy types will fail and therefore be ignored)


## [v0.1.1] - 05/04/2023

First release of the FLEPseq2 pipeline.

### Added
 - This changelog file to trace the different versions.
