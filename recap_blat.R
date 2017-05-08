#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
#
# Date: 170307
# Author: Gabriele Girelli
# Description:	converts a multi-sequence PSL file into a more readable format,
# 				with one sequence per row and the following columns:
# 				i: sequence ID number
# 				qName: sequence ID name (from PSL qName)
# 				k: sequence length in nt
# 				seq: sequence
# 				hX: number of genomic loci with X homology
# 					X from 1 to k
# 				maxHomology: highest homology with a genomic locus
#
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

library(argparser)
suppressMessages(library(data.table))
library(parallel)
library(readr)

# PARAMS =======================================================================

parser = arg_parser('Recap BLAT output.', name = 'recap_blat.R')

parser = add_argument(parser, arg = 'seq_dir',
	help = 'Directory containing the sequence file.',)
parser = add_argument(parser, arg = 'psl_dir',
	help = 'Directory containing the psl file.',)
parser = add_argument(parser, arg = 'out_dir',
	help = 'Output directory.',)

parser = add_argument(parser, arg = '--seq',
	help = 'Sequence file name.',
	default = 'bc.noHeader.fa')
parser = add_argument(parser, arg = '--psl',
	help = 'PSL file name.',
	default = 'merged_psl.psl')
parser = add_argument(parser, arg = '--out',
	help = "Output file name.",
	default = "merged_table.tsv")

parser = add_argument(parser, arg = '--nst',
	help = 'Number of sub tables.',
	default = 1000)
parser = add_argument(parser, arg = '--ncores',
	help = 'Number of cores for parallelization.',
	default = 30)

parser = add_argument(parser, arg = '--suff',
	help = "Output suffix.", default = '')


# Parse arguments --------------------------------------------------------------
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# PSL colnames based on PSL format definition
psl_colnames = c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert',
	'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize',
	'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount',
	'blockSizes', 'qStarts', 'tStarts')

if ("." != substr(suff, 1, 1) & 0 != nchar(suff)) suff = paste0(".", suff)

# RUN ==========================================================================

# Read sequences
cat(paste0(' · Reading sequence file...\n'))
seqt = read_delim(paste0(seq_dir, seq), '\t', col_names = c('seq'))$seq

# Read merged PSL table, and split it per query. Then, convert every sub-table
# in single rows and output them in single files. Finally, merge the sub-tables.
# REMEMBER to remove matches with inserts.
cat(paste0(' · Reading PSL file...\n'))
ccont = read_delim(paste0(psl_dir, psl), '\t', col_names = psl_colnames)
ccont = by(ccont, ccont$qName,
	FUN = function(x) {
		cat('Generated ', x$qName[1], ' sub-table.\n')
		return(x[x$tNumInsert == 0,])
	}
)

# A trick to use less memory: divide the sub-tables in groups of X sub-tables.
ccont = ccont
ccont2 = list()
n_tables = length(ccont)
start_group = seq(1, n_tables, nst)
i = 0
cat(paste0(' · Generating groups...\n'))
while(0 != length(ccont)) {
	i = i + 1
	start = start_group[i]
	end = min(start + nst - 1, n_tables)
	ccont2[[i]] = ccont[seq(1, nst)]
	ccont2[[i]][which(is.na(names(ccont2[[i]])))] = NULL
	ccont[seq(1, nst)] = NULL
	cat(paste0('   ', i, '/', length(start_group), '\n'))
}

# Iterate on the grouped sub-tables
cat(paste0(' · Recapping...\n'))
ktab = rbindlist(lapply(ccont2,
	FUN = function(ccont) {

		# Parallelize on the sub-tables
		ktab = rbindlist(mclapply(ccont,
			FUN = function(sub_ccont) {
				qName = sub_ccont$qName[1]
				cat(paste0(qName, '\n'))

				i = as.numeric(unlist(strsplit(qName, '_', fixed = T))[2])
				k = unique(sub_ccont$qSize)
				seq = seqt[i]
				out = data.frame(i, qName, k, seq)

				match_possible_values = seq(k)
				maxHomology = 0
				for (x in match_possible_values) {
					n = length(which(sub_ccont$matches == x))
					if (0 != n) maxHomology = max(maxHomology, x)
					out[,paste0('X', x)] = n
				}
				out$maxHomology = maxHomology

				return(out)
			}, mc.cores = ncores
		))
		cat(' · Writing output...\n')
		write.table(ktab, paste0(out_dir, out, suff), quote = F,
			row.names = F, col.names = F, sep = "\t", append = T)

		return(ktab)
	}
))

cat(' ~ END ~ ')

# END ==========================================================================

################################################################################
