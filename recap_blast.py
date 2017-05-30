#!/usr/bin/python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.0.1
# Date: 170530
# Description: Recap BLAST output in a more efficient format.
# Notes: BLAST output is relatively large (1e2 GB), thus go line by line.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(
	description = 'Recap BLAST output.'
)

# Add params
parser.add_argument('inFile', type = str, nargs = 1,
	help = 'Path to BLAST output with --outfmt 6.')
parser.add_argument('outFile', type = str, nargs = 1,
	help = 'Path to output recap file (will be created).')
parser.add_argument('-l', '--length', type = int, nargs = 1,
	help = 'Length in nt of the BLASTed oligomers.',
	default = [20])
parser.add_argument('-d', '--delimiter', type = str, nargs = 1,
	help = 'BLAST table delimiter.',
	default = ["\t"])

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
inFile = args.inFile[0]
outFile = args.outFile[0]
length = args.length[0]
sep = args.delimiter[0]

# Log to screen the settings
print("""
Input file: %s
Output file: %s
Oligomer length (nt): %i
Delimiter: '%s'
""" % (inFile, outFile, length, sep))

# FUNCTIONS ====================================================================

# RUN ==========================================================================

# Open pointer to outFile
fout = open(outFile, 'w')

# Print header
header = ['seq']
header.extend(["h" + str(i + 1) for i in range(length)])
header.append("hmax")
header = "\t".join(header) + "\n"
fout.write(header)

# Define default counter variables
cur_mer = 0
cur_counts = [0 for i in range(length)]
cur_max = 0

with open(inFile) as f:
	for line in f:
		row = line.strip().split(sep)

		# Check if the previous oligomeris concluded
		if not row[0] == cur_mer:
			if not cur_mer == 0:
				print(cur_mer)
				# Check max homology
				cur_max = max([i for i in range(length) if not cur_counts[i] == 0]) + 1

				# Prepare output
				outlist = [cur_mer]
				outlist.extend(cur_counts)
				outlist.append(cur_max)

				# Append recap to outFile if this is not the first oligo
				outstring = "\t".join([str(e) for e in outlist]) + "\n"
				fout.write(outstring)
			
			# Reset counts
			cur_mer = row[0]
			cur_counts = [0 for i in range(length)]
			cur_max = 0

			# Start new oligomer with default counts

		# Skip if it contains gaps
		if not 0 == int(row[5]):
			continue

		# Update oligomer counts otherwise
		# Read line, subtract mismatches from match length
		score = int(row[3]) - int(row[4])
		if score > length:
			print("!!!ERROR!!! Score [%i] is too big.\n%s" %(score, line))
		cur_counts[score - 1] += 1

# Close outFile pointer
fout.close()

# END ==========================================================================

################################################################################
