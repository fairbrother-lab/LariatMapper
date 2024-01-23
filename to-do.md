### Current top priority
- Change 3' end mapping to mapping to transcriptome followed by filtering by alignment location

### Definitely, sometime
- Move reference files + their args to one reference files folder and add a run-once script to prepare it
- Add unit tests
- Make this into a software package for installing as a command line command and distributing via anaconda/pip/etc

### Probably, sometime
- Add keep-intermediates optional argument to keep intermediate file
- Add individual arg-checks to larmap_rush.sh
- Move filter_lariats.py from larmap_run.sh into map_lariats.sh

### Maybe
- Find a way to map reads to 5'ss sequences instead of converting reads to an index and mapping 5'ss sequences to it (this is one of the most time-consuming steps in the pipeline)
- Change filter_threep_alignments.py to output all valid lariat mappings for each read instead of choosing 1 semi-randomly
- Make pipeline work with strand-specific sequencing data 
- Make pipeline work with 3rd-gen sequencing technology like Nanodrop/PacBio
- Refactor filter python scripts by replacing dicts and lists with pandas DataFrames (could be easier to read and maintain?)