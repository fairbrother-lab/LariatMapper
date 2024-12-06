suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(gt))



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
base_table_path = paste0(plots_dir, "Branchpoint_base_composition.html")
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

# Format values in table for gt
# Output will be <COUNT> \n <PERCENT>
cell_val = function(x) {
	perc = sprintf("%.0f%%", x / nrow(lariats) * 100)
	paste(x, perc, sep='<br>')
}
for (col in colnames(df)) {
	df[col] = as.character(apply(df[col], 2, cell_val))
}

# Add "In Read" column to act as the gt table row group label
df$type = "In Read"

# Create gt table
base_table <- gt(df, rownames_to_stub=TRUE, groupname_col = "type") %>%
	tab_header(md("Branchpoint base composition")) %>%
	# Spanner and column labels
	tab_spanner(label = md("In Genome"), columns = everything(),) %>%
	tab_style(style = cell_text(align="center", weight="bold", size=px(18)), locations = list(cells_column_spanners())) %>%
	tab_style(style = cell_text(align="center", weight="bold"), locations = list(cells_column_labels())) %>%
	# Stub and row labels
	tab_options(row_group.as_column = TRUE) %>%
	tab_style(style = cell_fill("#004D80"), locations = list(cells_stub(), cells_row_groups())) %>%
	tab_style(style = cell_text(align="center", v_align="middle", color="white", weight="bold"), locations = list(cells_stub())) %>%
    tab_style(style = cell_text(align="center", v_align="middle", color="white", weight="bold", size=px(18)),locations = list(cells_row_groups())) %>%
	# Fix row borders to match columns
    tab_style(style = cell_borders(sides="top", style="hidden"), locations = list(cells_stub(), cells_row_groups())) %>%
    tab_style(style = cell_borders(sides="right", style="solid", color='#89D3FE'), locations = list(cells_stub())) %>%
	# Whiten cells
	tab_style(style=cell_fill(color = "White"), locations = cells_body()) %>%
	# Table-wide formatting
	cols_width(
		c(A, C, G, T) ~ px(60),
		stub() ~ px(30),
		row_group() ~ px(90)
	) %>%
  	opt_stylize(style=2) %>%
	fmt_markdown() %>%
	tab_options(
		heading.title.font.size=px(24),
		heading.background.color="#003760",
		data_row.padding=px(12)
	) 
# Save plot to file
suppressMessages(gtsave(base_table, base_table_path))


# Density plot of branchpoint distance to 3'ss
# Create base plot
df = lariats[abs(lariats$bp_dist_to_threep)<=DIST_LIMIT,]
dists = (
	ggplot(df, aes(x=bp_dist_to_threep)) 
		+ geom_density(alpha=0.5, adjust=0.75) 
		+ scale_x_continuous(limits=c(-DIST_LIMIT,0), breaks=seq(-DIST_LIMIT, 0, 10))
		+ scale_y_continuous(labels=scales::percent)
		+ labs(x="Branchpoint distance to 3'ss (nt)", 
				y='Density')
		+ theme(plot.title = element_text(hjust=0.5, size=24),
				axis.text=element_text(size=14),
				axis.title=element_text(size=22),
		)
)
# Get peak height for placing bracket relative to peak
peak_height = max(ggplot_build(dists)$data[[1]]$y)
within_70_height = peak_height * 1.1
# within_limit_height = peak_height * 1.2
distal_height = peak_height * 1.25

dist_70 = sum(abs(lariats$bp_dist_to_threep)<=70)
dist_70_label = sprintf("%s lariat reads (%.0f%%) within 70nt", dist_70, dist_70 / nrow(lariats) * 100)
distal = sum(abs(lariats$bp_dist_to_threep)>DIST_LIMIT)
distal_label = sprintf("%s lariat reads (%.0f%%) beyond %snt", distal, distal / nrow(lariats) * 100, DIST_LIMIT)

# Add brackets
dists = (
	dists 
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