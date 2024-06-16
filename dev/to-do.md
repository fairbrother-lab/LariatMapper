### Current priorities
- Make it work and make sure its output is better than the original version

### Next priorities
- Find a way to map reads to 5'ss sequences instead of converting reads to an index and mapping 5'ss sequences to it (this is one of the most time-consuming steps in the pipeline)

### In the future
- Add version specifications to the dependencies in environment.yaml
- Convert map_lariats.sh into a python script for consistency and greater flexibility
- Option to input previously-created BAM as mapped_reads.bam, skipping the initial mapping step. Maybe make this the default?
- Change build_references.py so the copied reference files (fasta and repeatmasker) are always gzipped
- Add in the option to build_references.py to have it build the hisat2 index itself, maybe with a warning that it'll take a long time to do so
- Make this into a software package for installing as a command-line command and distributing via pip/pypi and/or anaconda
- Refine the repetitive region identification and filtering in filter_lariats.py, we should only consider cases where RNA transcribed from the region could be mistaken for a lariat due to it fitting the gapped inverted alignment and NOT a linear alignment
- Add support for stranded RNA-seq
- Add support for long-read sequencing (Nanopore/PacBio)

### Maybe
- Add a demo folder with scripts and reference files
- Create a test collection of reads corresponding to all possible read classification outputs (lariats, circularized introns, template-switching from lariat, normal mRNA/pre-MRNA from repetitive region, pre-mRNA, mRNA, lncRNA, intronic and starting at 5'ss, and jumbled nonsense sequence). Get the read seqs from the highest-confidence true-positive examples from literature (e.g. most commonly observed pre-mRNA species, circularized introns confirmed by sequencing and direct isolation, etc)
- Package this into a Docker image for distribution
- Add detection of reads from backspliced circular RNA, maybe using previously-published software as our implementation or as a validation method for an in-house implementation
