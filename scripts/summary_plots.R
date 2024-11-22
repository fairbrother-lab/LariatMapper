suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))



### Constants
BASES = c("A", "C", "G", "T")
BASE_COLORS = c("#008000","#0000ff", "#ffa600", "#ff0000")
DIST_LIMIT = 100
PLOT_DPI = 300



### Input
# Set args
args = commandArgs(trailingOnly = TRUE)
option_list = list(
    make_option(c("-o", "--output_base"), type = "character", default = NULL, 
                help = "Path to output", metavar = "OUT_PATH")
)

# Parse args
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

output_base <- opts$output_base

# Variable paths
lariats_file = paste(output_base, "lariat_reads.tsv", sep='')
plots_dir = paste0(output_base, "plots/")
pie_path = paste0(plots_dir, "Branchpoint_base_composition.png")
dists_path = paste0(plots_dir, "Branchpoint_threep_distance.png")

# Setup
theme_set(theme_classic())
if (!dir.exists(plots_dir)){
	dir.create(plots_dir)
}

# Check if lariats file exists
# If not, just exit
if (!file.exists(lariats_file)){
	quit(status=0)
}

# Load lariat datatable
lariats = read.table(lariats_file, sep='\t', header=TRUE)
lariats$read_bp_nt = factor(lariats$read_bp_nt, levels=BASES, ordered=TRUE)
lariats$genomic_bp_nt = factor(lariats$genomic_bp_nt, levels=BASES, ordered=TRUE)

# Check if there are any lariats
# If not, just exit
if (nrow(lariats) == 0){
	quit(status=0)
}



### Plots
# Pie chart of branchpoint base composition
# Calculate the frequencies for read_bp_nt
read_bp_freqs = as.data.frame(table(lariats$read_bp_nt))
read_bp_freqs$prop = read_bp_freqs$Freq / sum(read_bp_freqs$Freq)
read_bp_freqs$group = "In Read"

# Calculate the frequencies for genomic_bp_nt
genomic_bp_freqs = as.data.frame(table(lariats$genomic_bp_nt))
genomic_bp_freqs$prop = genomic_bp_freqs$Freq / sum(genomic_bp_freqs$Freq)
genomic_bp_freqs$group = "In Genome"

# Combine the two data frames
nt_freqs = rbind(read_bp_freqs, genomic_bp_freqs)
colnames(nt_freqs) = c("base", "freq", "prop", "group")
nt_freqs = nt_freqs[nt_freqs$freq>0,]
nt_freqs$base = factor(nt_freqs$base, levels=BASES, ordered=TRUE)
nt_freqs$group = factor(nt_freqs$group, levels=c("In Read", "In Genome"), ordered=TRUE)
nt_freqs$label = sprintf("%.0f%%", nt_freqs$prop * 100)
# nt_freqs$label = format(nt_freqs$freq, big.mark=',')
# nt_freqs$label = paste0(format(nt_freqs$freq, big.mark=','), " (", sprintf("%.0f%%", nt_freqs$prop * 100), ")")

pie = (
	ggplot(nt_freqs, aes(x="", y=freq, fill=base))
		+ geom_col(color="black")
		+ geom_text(aes(x=1.2, label=label), position=position_stack(vjust=0.5), size=6, color="White", fontface="bold") 
		+ coord_polar(theta='y') 
		+ scale_x_discrete(breaks=NULL) 
		+ scale_y_continuous(breaks=NULL) 
		+ scale_fill_manual(values=BASE_COLORS, guide=guide_legend(position='bottom')) 
		+ facet_wrap(~group) 
		+ labs(title="Branchpoint base composition", x="", y="", fill="") 
		+ theme(plot.title=element_text(hjust=0.5, size=24),
				strip.text=element_text(size=18),
				legend.text=element_text(size=18),
				panel.border=element_rect(color="black", fill=NA)
				)
)
suppressMessages(ggsave(pie_path, dpi=PLOT_DPI))


# Density plot of branchpoint distance to 3'ss
# Create base plot
df = lariats[abs(lariats$bp_dist_to_threep)<=DIST_LIMIT,]
dists = (
	ggplot(df, aes(x=bp_dist_to_threep)) 
		+ geom_density(alpha=0.5, adjust=0.75) 
		+ scale_x_continuous(limits=c(-DIST_LIMIT,0), breaks=seq(-DIST_LIMIT, 0, 10))
		+ labs(x="Branchpoint distance to 3'ss (nt)", 
				y='Density')
		+ theme(plot.title = element_text(hjust=0.5, size=24),
				axis.text=element_text(size=14),
				axis.title=element_text(size=22),
		)
)
# Get peak height for placing bracket relative to peak
peak_height = max(ggplot_build(dists)$data[[1]]$y)
dist_70_height = peak_height * 1.1
dist_limit_height = peak_height * 1.2

dist_70_prop = sum(abs(lariats$bp_dist_to_threep)<=70)/nrow(lariats)
dist_70_label = paste0(sprintf("%.1f%%", dist_70_prop * 100), ' of lariat reads')
dist_limit_prop = sum(abs(lariats$bp_dist_to_threep)<=DIST_LIMIT)/nrow(lariats)
dist_limit_label = paste0(sprintf("%.1f%%", dist_limit_prop * 100), ' of lariat reads')

# Add brackets
dists = (
	dists 
	# dist 70 bracket
		+ annotate("segment", x = -70, xend = 0, y = dist_70_height, yend = dist_70_height) 
		+ annotate("segment", x = -70, xend = -70, y = dist_70_height*0.98, yend = dist_70_height) 
		+ annotate("segment", x = 0, xend = 0, y = dist_70_height*0.98, yend = dist_70_height)
		+ annotate("text", x = -35, y = dist_70_height*1.02, label = dist_70_label, size = 5, hjust=0.5, vjust=0) 
	# dist limit bracket
		+ annotate("segment", x = -DIST_LIMIT, xend = 0, y = dist_limit_height, yend = dist_limit_height) 
		+ annotate("segment", x = -DIST_LIMIT, xend = -DIST_LIMIT, y = dist_limit_height*0.98, yend = dist_limit_height) 
		+ annotate("segment", x = 0, xend = 0, y = dist_limit_height*0.98, yend = dist_limit_height)
		+ annotate("text", x = -DIST_LIMIT/2, y = dist_limit_height*1.02, label = dist_limit_label, size = 5, hjust=0.5, vjust=0)
)
suppressMessages(ggsave(dists_path, dpi=PLOT_DPI))