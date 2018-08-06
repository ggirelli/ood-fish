#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
#
# Date: 170307
# Author: Gabriele Girelli
# Description:	filter HDFE table and produce PMC compatible format.
#
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

library(argparser)
suppressMessages(library(data.table))
library(readr)

# PARAMS =======================================================================

parser = arg_parser("Filter BLAT recap table based on SDFE and maxHomology.",
	name = "hdfe_filter.R")

parser = add_argument(parser, arg = "in_dir",
	help = "Input directory.")
parser = add_argument(parser, arg = "out_dir",
	help = "Output directory.")

parser = add_argument(parser, arg = "--fin",
	help = "Input file name.",
	default = "hdfe.dat")
parser = add_argument(parser, arg = "--out",
	help = "Output file name.",
	default = "hdfe_filtered.dat")
parser = add_argument(parser, arg = "--pmc",
	help = "Output PMC file name.",
	default = "hdfe_filtered_graph.txt")

parser = add_argument(parser, arg = "--ht",
	help = "Homology threshold (positive).",
	default = -9)

# Parse arguments --------------------------------------------------------------
p = parse_args(parser)

# Attach argument values to variables
attach(p["" != names(p)])

# Updates
ht = -abs(ht)

# Print settings
cat(paste0("HDFE threshold: ", ht, " kcal/mol\n"))

# RUN ==========================================================================

# Read input.
cat(" · Reading input...\n")
t = read_delim(paste0(in_dir, fin), ' ',
	col_names = c('start', 'end', 'hdfe'), col_types = 'iid')

# Filter.
cat(" · Filtering...\n")
tf = t[t$hdfe >= ht,]

# Write filtered table.
cat(" · Writing table...\n")
write.table(tf, paste0(out_dir, out),
	quote = F, col.names = F, row.names = F, sep = ' ')

# Prepare PMC-compatible table.
g = tf[, c(1, 2)]

# Write PMC output.
cat(" · Writing graph...\n")
write.table(g, paste0(out_dir, pmc),
	quote = F, col.names = F, row.names = F, sep = ' ')

cat(" ~ END ~ ")

# END ==========================================================================

################################################################################