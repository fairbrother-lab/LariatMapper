| read name | details |
|--|--|
| mRNA	| chr1, DDX11L1, exons #3-#4, head is 638-688 and tail is 965-1015, BP at 49 in read |
| pre-mRNA	| chr1, DDX11L1, intron #4 to exon #5 junction, 1200-1300 |
| true_lariat_1	| chr1, DDX11L1, intron #3, head is 892-942 and tail is 688-738, BP at 49 in read |
| true_lariat_2	| chr1, DDX11L1, intron #5, head is 1374-1424 and tail is 1365-1415, BP at 49 in read |
| fivep_fail_5bp_up_mismatch | derived from true_lariat_1, only kept 5'ss and added the 5nt upstream in genome |
| fivep_fail_enough_head_seq | derived from true_lariat_1, only kept 5'ss and moved it to 5nt from the read's start |  
| head_fail_mismatches | derived from true_lariat_1, change the last 6 A's in head to T's |
| head_fail_gaps | derived from true_lariat_1, delete TTTCCC from head and append the next 6nt to its end |
| temp_switch_1 | derived from true_lariat_1, move head to upstream of GTGAG in chr1 |  
| temp_switch_2 | derived from true_lariat_2, move head to upstream of GTTCA in chr2 |
| head_fail_overlap_introns | derived from true_lariat_1, move to start of exon #5 |
| head_fail_fivep_intron_match | derived from true_lariat_1, move head to 150-200 in chr2 |
| head_fail_fivep_intron_match | derived from true_lariat_1, reverse-complement head |
| head_fail_fivep_intron_match | derived from true_lariat_1, move head to 700-750 |
| circular | derived from true_lariat_1, replace head with seq at end of intron #3 |
| repeat | chr2, FOXN2, intron #1, head 200-250 and tail 133-183 |
Positions are 0-based inclusive start, exclusive end