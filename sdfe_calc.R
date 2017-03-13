#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
#
# Date: 170307
# Author: Gabriele Girelli
# Description:	calculate self dimerization free energy
# 				and add it to a BLAT recap table.
#
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

library(argparser)
suppressMessages(library(data.table))
library(parallel)
library(readr)

source('/home/gire/Desktop/BiCro/Code/repo/00_gabriele/Sequencing/seq_lib.R')

# PARAMS =======================================================================

parser = arg_parser('Calculate self dimerization free energy.',
	name = 'sdfe_calc.R')

parser = add_argument(parser, arg = 'rec_dir',
	help = 'Directory containing the recap table.',)
parser = add_argument(parser, arg = 'out_dir',
	help = 'Output directory.',)

parser = add_argument(parser, arg = '--rec',
	help = 'Recap table file name.',
	default = 'merged_table.tsv')
parser = add_argument(parser, arg = '--out',
	help = "Output file name.",
	default = "merged_table.sdfe.tsv")

parser = add_argument(parser, arg = '--ncores',
	help = 'Number of cores for parallelization.',
	default = 30)

parser = add_argument(parser, arg = '--suff',
	help = "Output suffix.", default = '')

# Parse arguments --------------------------------------------------------------
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# RUN ==========================================================================

# Calculate self dimerization energy for every query in the 
cat(paste0(' · Reading recap table...\n'))
t = read_delim(paste0(rec_dir, rec), '\t')

# Initialize dimerization suite and start calculation
ds = dimerization_suite(ncores = 1)

# Calculate SDFE
cat(paste0(' · Calculating self-dimerization free energy of ', nrow(t), ' ',
	t$k[1], '-mers... [using ', ncores, ' cores]\n'))
l = unlist(mclapply(t$seq, FUN = ds$get_self_dimer_free_energy,
	only.top = T, mc.cores = ncores))

# Add SDFE column
t$sdfe = l

# Write output
cat(paste0(' · Writing output...\n'))
write.table(t, paste0(out_dir, out, suff),
	quote = F, row.names = F, col.names = T, sep = '\t')

cat(' ~ END ~ ')

# END ==========================================================================

################################################################################