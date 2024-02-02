### Current top priority
- Change 3' end mapping to mapping to transcriptome followed by filtering by alignment location

### Definitely, sometime
- Move reference files + their args to one reference files folder and add a run-once script to prepare it (replacing get_splice_site_seqs.py)
- Make this into a software package for installing as a command-line command and distributing via anaconda/pip/etc
- Update larmap_setup.py to work with current version or remove it
- Update larmap_merge.py to work with current version or remove it
- Make output base name into option, default to no prefix

### Probably, sometime
- Add keep-intermediates optional argument to keep intermediate files
- Move filter_lariats.py from larmap_run.sh into map_lariats.sh
- Change python scripts print statements to logs and add logging calls that only print in debug mode (e.g. input arg values, DataFrame heads)
- Add unit tests (for each filter, at least 1 artificial read that reaches it and fails it)
- Specify package versions in environment.yaml to future-proof 
- Add individual arg-checks to larmap_rush.sh
- Refactor filter_fivep_alignments.py to replace dicts with pandas DataFrames

### Maybe
- Find a way to map reads to 5'ss sequences instead of converting reads to an index and mapping 5'ss sequences to it (this is one of the most time-consuming steps in the pipeline)
- Change filter_threep_alignments.py to output all valid lariat mappings for each read instead of choosing 1 semi-randomly
- Make pipeline work with strand-specific sequencing data 
- Make pipeline work with 3rd-gen sequencing technology like Nanodrop/PacBio
- Refactor filter python scripts by replacing dicts and lists with pandas DataFrames (could be easier to read and maintain?)