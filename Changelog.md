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

First stable release of the pipeline.

## [v0.2.0] 15/04/2023

### Added
 - a log for the merge_read_info.R to track the different steps (and the number of lines of each step)
 - additional tail detection even when there is no polyA tail (some transcripts like the elongating_3_mapping_low_accuracy types will fail and therefore be ignored)


## [v0.1.1] - 05/04/2023

First release of the FLEPseq2 pipeline.

### Added
 - This changelog file to trace the different versions.
