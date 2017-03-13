#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
#
# Date: 170307
# Author: Gabriele Girelli
# Description:	filter a BLAT recap table based on maxHomology and sdfe.
#
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

library(argparser)
suppressMessages(library(data.table))
library(readr)

# PARAMS =======================================================================

parser = arg_parser("Filter BLAT recap table based on SDFE and maxHomology.",
	name = "sdfe_filter.R")

parser = add_argument(parser, arg = "in_dir",
	help = "Input directory.")
parser = add_argument(parser, arg = "out_dir",
	help = "Output directory.")

parser = add_argument(parser, arg = "--fin",
	help = "Input file name.",
	default = "merged_table.sdfe.tsv")
parser = add_argument(parser, arg = "--out",
	help = "Output file name.",
	default = "merged_table.filtered.sdfe.tsv")

parser = add_argument(parser, arg = "--ht", nargs = 1,
	help = "Homology threshold.",
	default = .75, type = class(.0))
parser = add_argument(parser, arg = "--st", nargs = 1,
	help = "Self-dimerization free energy threshold (positive).",
	default = -5, type = class(0))

# Parse arguments --------------------------------------------------------------
p = parse_args(parser)

# Attach argument values to variables
attach(p["" != names(p)])

# Updates
st = -abs(st)

# Print settings
cat(paste0("Homology threshold: ", ht, "%\n"))
cat(paste0("SDFE threshold: ", st, " kcal/mol\n"))

# FUNCTIONS ====================================================================

add_trailing_slash = function(s, del = "/") {
	ls = nchar(s)
	if ( del != substr(s, ls, ls) ) s = paste0(s, del)
	return(s)
}

# RUN ==========================================================================

# Add trailing slashes
in_dir = add_trailing_slash(in_dir)
out_dir = add_trailing_slash(out_dir)

# Read input table
cat(" · Reading recap table with SDFE...\n")
t = read_delim(paste0(in_dir, fin), "\t")

# Calculate threshold
thr = t$k[1] * ht

# Apply threshold
cat(" · Applying threshold...\n")
tf = t[t$maxHomology <= thr & t$sdfe > st,]

# Output
cat(" · Writing output...\n")
write.table(tf, paste0(out_dir, out),
	quote = F, row.names = F, col.names = T, sep = "\t")

cat(" ~ END ~ \n")

# END ==========================================================================

################################################################################