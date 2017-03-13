#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
#
# Date: 170307
# Author: Gabriele Girelli
# Description:	merge split BLAT recap tables.
#
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

library(argparser)
suppressMessages(library(data.table))
library(readr)

# PARAMS =======================================================================

parser = arg_parser('Merge BLAT recaps.', name = 'recap_merge.R')

parser = add_argument(parser, arg = 'recap_dir',
	help = 'Directory containing the recap files.',)
parser = add_argument(parser, arg = 'out_dir',
	help = 'Output directory.',)

parser = add_argument(parser, arg = '--out',
	help = "Output file name.",
	default = "merged_table.tsv")

# Parse arguments --------------------------------------------------------------
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# RUN ==========================================================================

# Identify split recaps
recaps = list.files(recap_dir)

# Read and merge split recaps
cat(' · Read and merge split recap tables.\n')
recaps = rbindlist(lapply(paste0(recap_dir, recaps), FUN = read_delim, '\t',
	col_names = F))
colnames(recaps) = c('i', 'qName', 'k', 'seq',
	paste0('h', seq(recaps$X3[1])), 'maxHomology')

# Fix NA IDs
cat(' · Fixing eventual NA IDs.\n')
nais = which(is.na(recaps$i))
recaps$i[nais] = unlist(lapply(nais,
	FUN = function(i) {
		as.numeric(unlist(strsplit(recaps$qName[i], '_', fixed = T))[2])
	}
))

# Fix duplications caused by merging
cat(' · Fixing duplications caused by merge.\n')
dup_qNames = recaps$qName[which(duplicated(recaps$qName))]
dedup = rbindlist(lapply(dup_qNames,
	FUN = function(dname) {
		subt = as.data.frame(recaps[recaps$qName == dname,])
		new_row = subt[1,]
		new_row[, paste0('h', seq(new_row$k))] = apply(
			subt[, paste0('h', seq(new_row$k))], FUN = sum, MARGIN = 2)
		new_row$maxHomology = max(subt$maxHomology)

		return(new_row)
	}
))
recaps = recaps[-which(recaps$qName %in% dup_qNames),]
recaps = rbind(recaps, dedup)
recaps = recaps[order(recaps$i),]

# Output
cat(' · Writing output...\n')
write.table(recaps, paste0(out_dir, out),
	quote = F, row.names = F, col.names = T, sep = '\t')

cat(' ~ END ~ ')

# END ==========================================================================

################################################################################
