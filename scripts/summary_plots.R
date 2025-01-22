suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(scales))



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
base_table_path = paste0(plots_dir, "Branchpoint_base_composition.png")
dists_path = paste0(plots_dir, "Branchpoint_threep_distance.png")

# Setup
theme_set(theme_classic())
if (!dir.exists(plots_dir)){
	dir.create(plots_dir)
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



# Derive base composition table 
df = table(lariats[, c("read_bp_nt", "genomic_bp_nt")]) 
df = as.data.frame.matrix(df)

# Format values in table for ggplot
df <- df %>%
  rownames_to_column(var = "In_Read") %>%
  pivot_longer(cols = -In_Read, names_to = "In_Genome", values_to = "Count") %>%
  mutate(Percent = Count / nrow(lariats),
  		In_Read = factor(In_Read, levels = rev(BASES), ordered = T),
		In_Genome = factor(In_Genome, levels = BASES, ordered = T)
)

# Create table like a heatmap
base_table <- ggplot(df, aes(x = In_Genome, y = In_Read, fill = Percent)) +
  geom_tile(color="#003760", linewidth=1.5) +
  geom_text(aes(label = Count), color = "black", size=6) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "white", high = "#4daf4a", labels=percent_format()) +
  labs(x = "In Genome", y = "In Read") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 30, face = "bold", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
	legend.title = element_blank(),
	legend.text = element_text(size = 18, color = "white"),
	legend.key.size = unit(30, 'pt'),
	legend.key.height = unit(50, 'pt'),
    plot.background = element_rect(fill = "#003760"),
	plot.margin = margin(40, 40, 40, 40, "pt"),
    axis.title = element_text(size = 48, face = "bold", color = "white"),
  )

# Save plot to file
suppressMessages(ggsave(base_table_path, base_table, height=8, width=8))



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
densities = ggplot_build(dists)$data[[1]]$y
peak_height = ifelse(is.null(densities) || any(is.na(densities)), 0.8, max(densities))
within_70_height = peak_height * 1.1
distal_height = peak_height * 1.25

dist_70 = sum(abs(lariats$bp_dist_to_threep)<=70)
dist_70_label = sprintf("%s lariat reads (%.0f%%) within 70nt", dist_70, dist_70 / nrow(lariats) * 100)
distal = sum(abs(lariats$bp_dist_to_threep)>DIST_LIMIT)
distal_label = sprintf("%s lariat reads (%.0f%%) beyond %snt", distal, distal / nrow(lariats) * 100, DIST_LIMIT)

# Add brackets
dists = (
	dists 
		+ scale_y_continuous(limits=c(0, distal_height), labels=scales::percent)
	# dist 70 bracket
		+ annotate("segment", x = -70, xend = 0, y = within_70_height, yend = within_70_height) 
		+ annotate("segment", x = -70, xend = -70, y = within_70_height*0.98, yend = within_70_height) 
		+ annotate("segment", x = 0, xend = 0, y = within_70_height*0.98, yend = within_70_height)
		+ annotate("text", x = -35, y = within_70_height*1.02, label = dist_70_label, size = 5, hjust=0.5, vjust=0) 
	# distal branchpoints arrow
		+ annotate("segment", x = -DIST_LIMIT*0.9, xend = -DIST_LIMIT, y = distal_height, yend = distal_height, arrow=arrow())
		+ annotate("text", x = -DIST_LIMIT*0.89, y = distal_height, label = distal_label, size = 5, hjust=0, vjust=0.5)
)
# Save plot to file
suppressMessages(ggsave(dists_path, dpi=PLOT_DPI))