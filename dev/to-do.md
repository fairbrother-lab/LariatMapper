### Current priority
- Find a way to map reads to 5'ss sequences instead of converting reads to an index and mapping 5'ss sequences to it (this is one of the most time-consuming steps in the pipeline)
- Make it work

### Next priority


### In the future
- Start using Releases for updates
- Add a demo folder with scripts and reference files
- Change build_references.py so the copied reference files (fasta and repeatmasker) are always gzipped 
- Add support for stranded RNA-seq
- Add support for long-read sequencing (Nanopore/Pacbio)
- Add version specifications to the dependencies in environment.yaml
- Option to input previously-created BAM as mapped_reads.bam, skipping the initial mapping step. Maybe make this the default?

### Maybe
- Add unit testing for development
- Make this into a software package for installing as a command-line command and distributing via pip/pypi and/or anaconda
- Package this into a Docker image
- Add detection of reads from backspliced RNA
- Instead of reading all alignments in filter_fivep_alignments.py and filter_trimmed_alignments.py into one object then passing chunks to processes, count the number of recorded alignments and assign a range to each process, letting each process read alignments and process them on their own. This would also allow processing smaller chunks in case there are too many alignments to have them all in memory at once.

