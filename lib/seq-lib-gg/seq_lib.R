# @description: Base library for sequencing data analysis.
# @author: Gabriele Girelli
# @version: 0.1.0



# DEPENDENCIES =================================================================

library(ggplot2)
library(parallel)
library(reshape2)
library(XML)

# FUNCTIONS ====================================================================

rev_str = function(s) {
	# Reverse a string
	
	paste(rev(unlist(strsplit(s, '', fixed=T))), collapse='')
}

complement = function(na, t='DNA') {
	# Provided a nucleic acid, returns the complement.
	# Mus specify in `t` if the provided na is DNA or RNA.
	
	# Identify type
	t <- tolower(t)

	# Select alphabet
	if ( 'dna' == t ) {
		ab <- c('ATCG', 'TAGC')
	} else if ( 'rna' == t ) {
		ab <- c('AUCG', 'UAGC')
	} else {
		cat('ERROR: unknown na type.\n')
		return()
	}

	rab <- unlist(strsplit(tolower(ab[2]), ''))
	ab <- unlist(strsplit(tolower(ab[1]), ''))

	# Check provided string
	na <- tolower(na)
	for ( c in unlist(strsplit(na, '')) ) {
		if ( !c %in% ab ) {
			print(c)
			print(ab)
			cat('ERROR: the provided string conflicts',
				'with the selected alphabet.\n')
			return()
		}
	}

	# Calculate complement
	compl <- c()
	for ( c in unlist(strsplit(na, '')) ) {
		compl <- append(rab[which(ab == c)], compl)
	}
	compl <- toupper(paste(compl, collapse=''))

	return(rev_str(compl))
}

rc = function(na, t='DNA') {
	# Provided a nucleic acid, returns the reverse complement.
	# Mus specify in `t` if the provided na is DNA or RNA.

	return(rev_str(complement(na, t)))
}

get_base_pdistr = function(seql,
	ncores = 1,
	col = c('#74C957', '#0068CB', '#CC2343', '#009688', '#FF9800')
) {
	# Calculates the frequency of each base at each position from a list of
	# sequences of same length
	# 
	# Depends on: ggplot2, parallel, reshape2

	# Check that the sequences in the provided list have the same length
	seq_lens <- unique(unlist(lapply(seql, nchar)))
	if ( 1 != length(seq_lens) ) {
		cat('ERROR: found sequences of different length...\n')
		print(seq_lens)
		return(NULL)
	}

	# Make base matrix
	mb <- do.call(rbind, mclapply(seql,
		FUN=function(x) {
			unlist(strsplit(x, '', fixed = T))
		}
		, mc.cores = ncores
	))

	# Count and calculate frequencies
	mt <- apply(mb, MARGIN = 2, FUN = table)
	mt <- mt / length(seql)

	# From table to dataframe
	md <- data.frame(t(mt))
	md$position <- 1:nrow(md)

	# Plot
	p <- ggplot(
		data = melt(md, id.vars = 'position'),
		aes(x = position, y = value, colour = variable)
	) + geom_line() + ylim(c(0, 1))
	p <- p + labs(
		x = 'base position [nt]',
		y = 'base probability',
		title = 'Position-based base probability distribution'
	)
	p <- p + scale_colour_manual('',
		breaks = colnames(md)[-ncol(md)],
		labels = colnames(md)[-ncol(md)],
		values = col[1:(ncol(md) - 1)]
	)
	p <- p + geom_hline(
		yintercept = 0.25,
		color = '#989898',
		linetype = 'dashed'
	)
	print(p)

	#Output
	return(list(table = mt, plot = p))
}

# CLASSES ======================================================================

# IDT_tools --------------------------------------------------------------------

# Contains functions to query IDT OligoAnalyzer web services.
# Depends on: XML
# 
IDT_tools = function() {

	# Generate class
	it <- list(

		gsdfe = function(seq) {
			# Get Self-Dimerization Free Energy
			xmluri <- paste0('http://www.idtdna.com/AnalyzerService/',
				'AnalyzerService.asmx/SelfDimer?Sequence=', seq)
			xmlfile <- xmlTreeParse(xmluri)
			xmltop = xmlRoot(xmlfile)
			seqdata <- xmlSApply(xmltop, function(x) xmlSApply(x, xmlValue))
			min(as.numeric(seqdata$DimerInfo))
		},

		ghdfe = function(seq, lseq) {
			# Get Hetero-Dimerization Free Energy
			unlist(lapply(1:length(lseq),
				FUN=function(i) {
					cat(paste0(' · Hetero-dimerizing against seq #', i, '\n'))
					seq1 <- seq
					seq2 <- lseq[i]
					xmluri <- paste0('http://www.idtdna.com/AnalyzerService/',
						'AnalyzerService.asmx/HeteroDimer?Sequence=', seq1,
						'&SecondarySequence=', seq2)
					xmlfile <- xmlTreeParse(xmluri)
					xmltop <- xmlRoot(xmlfile)
					seqdata <- xmlSApply(xmltop,
						function(x) xmlSApply(x, xmlValue))
					min(as.numeric(seqdata$DimerInfo))
				}
			))
		},

		ghdfe_rc = function(seq, lseq) {
			# Get Self-Dimerization Free Energy with the reverse complement
			# of the provided sequences
			seq <- rc(seq, 'dna')
			unlist(lapply(1:length(lseq),
				FUN=function(i) {
					cat(paste0(' · Hetero-dimerizing against seq #', i, '\n'))
					seq1 <- seq
					seq2 <- lseq[i]
					xmluri <- paste0('http://www.idtdna.com/AnalyzerService/',
						'AnalyzerService.asmx/HeteroDimer?Sequence=', seq1,
						'&SecondarySequence=', seq2)
					xmlfile <- xmlTreeParse(xmluri)
					xmltop <- xmlRoot(xmlfile)
					seqdata <- xmlSApply(xmltop,
						function(x) xmlSApply(x, xmlValue))
					min(as.numeric(seqdata$DimerInfo))
				}
			))
		},

		gmt = function(seq,
			target_type = 'DNA',
			oligo_conc = 0.25,		# uM
			na_conc = 50,			# mM
			mg_conc = 0,			# mM
			dntp_conc = 0,			# mM
			verbose = F
		) {
			# Depends on XML
			xmluri <- paste0('http://www.idtdna.com/AnalyzerService/',
				'AnalyzerService.asmx/', 'Analyze?Sequence=', seq,
				'&TargetType=', target_type, '&OligoConc=', oligo_conc,
				'&NaConc=', na_conc, '&MgConc=', mg_conc,
				'&dNTPsConc=', dntp_conc)
			xmlfile <- xmlTreeParse(xmluri)
			xmltop = xmlRoot(xmlfile)
			Tm <- as.numeric(xmlSApply(xmltop,
				function(x) xmlSApply(x, xmlValue))$MeltTemp)
			if ( verbose )
				cat(paste0('IDT-calculated Tm: ', Tm, ' °C',
					' [', seq, ']\n'))
			return(Tm)
		}

	)

	# Assign class
	class(it) <- 'IDT_tools'

	# Output
	return(it)
}

# simple_match -----------------------------------------------------------------

# Provided two strings of equal length, return the match score
# by comparing the i-th char of the two strings.
# Depends on: parallel
# 
simple_match = function(s1 = '', s2 = '', ncores = 1) {

	# Declare default class
	m <- list(
		s1 = paste(s1, collapse=''),
		s2 = paste(s2, collapse=''),
		vmatch = NULL,
		score = 0,
		match_string = '',
		contig = list(
			score = 0,
			pos = NULL,
			seq = NULL
		)
	)
	class(m) <- 'simple_match'
	
	# Check that the strings have the same length
	if ( nchar(s1) != nchar(s2) ) {
		cat('[simple_match] ERROR: the provided strings',
			'must be of equal length.\n')
		m$score = -1
		return(m)
	}

	if ( 0 == nchar(s1) ) return(m)

	# Convert the sequences in vectors
	s1 <- unlist(strsplit(s1, '', fixed=T))
	s2 <- unlist(strsplit(s2, '', fixed=T))

	# Produce binary-match vector
	vmatch <- unlist(mclapply(1:length(s1),
		FUN=function(i) {
			if ( identical(s1[i], s2[i]) ) return(1)
			else return(0)
		}
		, mc.cores = ncores
	))

	# Match score
	score_simple <- sum(vmatch)

	# Biggest contiguous match
	sub_vmatch <- unlist(strsplit(paste(vmatch, collapse=''), '0', fixed=T))
	len_contig <- nchar(sub_vmatch)
	score_contig_max <- max(len_contig)
	len_contig <- len_contig[len_contig != 0]

	# Generate match string
	match_string <- rep(' ', length(s1))
	match_string[which( 1 == vmatch )] <- '·'

	if ( 0 == score_contig_max ) {
		pos_contig <- NULL
		seq_contig <- NULL
	} else {
		# Retrieve contiguous match positions
		pos_contig <- as.numeric(gregexpr(
			'01', paste(c(0, vmatch), collapse='')
		)[[1]])

		# Retrieve contiguous match sequences
		seq_contig <- unlist(mclapply(1:length(pos_contig),
			FUN=function(i) {
				pos <- pos_contig[i]
				len <- len_contig[i]
				paste(s1[pos:(pos + len - 1)], collapse='')
			}
			, mc.cores = ncores
		))

		# Update match string
		match_string[pos_contig[1]:(pos_contig[1] + len_contig[1] - 1)] <- '|'
	}

	# Generate class
	m <- list(
		s1 = paste(s1, collapse=''),
		s2 = paste(s2, collapse=''),
		vmatch = vmatch,
		score = score_simple,
		match_string = paste(match_string, collapse=''),
		contig = list(
			score = score_contig_max,
			len = len_contig,
			pos = pos_contig,
			seq = seq_contig
		)
	)
	class(m) <- 'simple_match'

	# Output
	return(m)
}

print.simple_match = function(sm) {
	# Prints a simple_match output in a readable way
	if ( !is.null(sm$contig$fe) ) {
		cat(paste0('\nScore: ', sm$score, '\nMax contig: ', sm$contig$score,
			'\nFree energy: ', min(sm$contig$fe, na.rm=T),
			'\n\n', sm$s1, '\n', sm$match_string, '\n', sm$s2, '\n'))
	} else {
		cat(paste0('\nScore: ', sm$score, '\nMax contig: ', sm$contig$score,
			'\n\n', sm$s1, '\n', sm$match_string, '\n', sm$s2, '\n'))
	}
}

print.simple_match_list = function(sml) {
	# Print the align_match output in a readable way
	if ( 0 == length(sml) ) return()
	l <- lapply(1:length(sml),
		FUN=function(i) {
			cat(paste0('\n\nMatch #', i))
			print(sml[[i]])
		}
	)
}

order.simple_match_list.score = function(ml) {
	order(unlist(lapply(ml, FUN=function(x) x$score)))
}

order.simple_match_list.contig = function(ml) {
	order(unlist(lapply(ml, FUN=function(x) x$contig$score)))
}

order.simple_match_list.free_energy = function(ml) {
	order(abs(unlist(lapply(ml, FUN=function(x) {
		if ( !is.null(x$contig$fe) ) {
			return(min(x$contig$fe, na.rm=T))
		} else {
			return(0)
		}
	}))))
}

order.simple_match_list = function(ml, by = 'score', decreasing = F) {
	# Order a simple_match_object.
	# Can be ordered by: score, contig, free_energy

	order <- switch(by,
		'score' = order.simple_match_list.score(ml),
		'contig' = order.simple_match_list.contig(ml),
		'free_energy' = order.simple_match_list.free_energy(ml)
	)

	if ( decreasing ) return(rev(order))
	return(order)
}

# align_fast -------------------------------------------------------------------

# Provided two strings, find the best alignments of the shortest one
# in the longest and return them alongside the score.
# Depends on: parallel
# 
align_fast = function(
	s1 = '', s2 = '',
	t = 'DNA',
	dg_table = NULL,
	self.dimer = F, hetero.dimer = F,
	ncores = 1
) {

	# Generate class
	af = list(
		s1 = s1, s2 = s2,				# Strings
		matches = c(),					# Sores the output
		t = t,							# Nucleic acid type (DNA|RNA)
		dg_table = dg_table,			# Standard dimerization dG table (path)
		self.dimer = self.dimer,		# Perform self-dimerization of s1
		hetero.dimer = hetero.dimer,	# Perform dimerization of s1 on s2
		ncores = ncores,				# Number of corse for parallel computing

		align = function() {
			# Revert if dimer
			if ( af$self.dimer ) af$s2 <- rc(af$s1, af$t)
			if ( af$hetero.dimer ) af$s2 <- complement(af$s2, af$t)

			# Identify long/short string
			s <- af$assign(af$s1, af$s2)

			# No null strings allowed
			if ( 0 == nchar(af$s1) | 0 == nchar(af$s2) ) return(af)

			# Run matching on begin/inner/end part
			af$matches <- c(
				af$begin(af$s1, af$s2),
				af$inner(af$s1, af$s2),
				af$end(af$s1, af$s2)
			)
			class(af$matches) <- 'simple_match_list'

			return(af)
		},

		assign = function(s1, s2) {
			# Order provided strings based on length
			if ( nchar(s1) >= nchar(s2) )
				return(list(long = s1, short = s2))
			else
				return(list(long = s2, short = s1))
		},

		revert_complement = function(m) {
			# Revert complement
			m$s2 <- rc(m$s2, af$t)
			m$contig$seq <- unlist(lapply(m$contig$seq, FUN=rc, af$t))
			return(m)
		},
		cal_contig_free_energy = function(m) {
			# Calculate dimerization free energy of every contig
			if ( 2 <=  m$contig$score ) {

				id_contig <- which(2 <= m$contig$len)
				fe_contig <- rep(NA, length(m$contig$len))

				fe_contig[id_contig] <- unlist(lapply(id_contig,
					FUN=function(i) {
						dimerization_suite()$gfe(m$contig$seq[i])
					}
				))

				m$contig$fe <- fe_contig

				# Re-generate match string
				match_string <- unlist(strsplit(m$match_string, '', fixed=T))
				match_string['|' == match_string] <- '·'
				top_cont <- which.min(fe_contig)
				pos <- m$contig$pos
				len <- m$contig$len
				chars_contig <- pos[top_cont]:(pos[top_cont] + len[top_cont] -1)
				match_string[chars_contig] <- '|'

				m$match_string <- paste(match_string, collapse='')
			} else {
				m$contig$fe <- 0	
			}

			return(m)
		},

		begin = function(s1, s2) {
			# Provided two strings, finds the best alignments of the shortest
			# one in the beginning of the longest and returns match score,
			# sequence and position
			
			# Identify long/short string
			s <- af$assign(s1, s2)
			len_short <- nchar(s$short)

			# Slide substrings of the short seq on the beginning of the long one
			lmatch <- mclapply(1:(len_short - 1),
				FUN=function(i) {
					sub_long <- substr(s$long, 1, i)
					sub_short <- substr(s$short, len_short - i + 1, len_short)

					m <- simple_match(sub_long, sub_short)
					m$pos_long <- c(1, i)
					m$pos_short <- c(len_short - i + 1, len_short)

					if ( self.dimer | hetero.dimer ) {
						m <- af$cal_contig_free_energy(m)
						m <- af$revert_complement(m)
					}

					return(m)
				}
				, mc.cores = af$ncores
			)

			return(lmatch)
		},
		inner = function(s1, s2) {
			# Provided two strings, finds the best alignments of the shortest
			# one in the longest and returns match score, sequence and position
			
			# Identify long/short string
			s <- af$assign(s1, s2)
			len_short <- nchar(s$short)
			len_long <- nchar(s$long)
			
			# Slide substrings of the short seq on the long one
			lmatch <- mclapply(1:(len_long - len_short + 1),
				FUN=function(i) {
					sub_long <- substr(s$long, i, i + len_short - 1)
					sub_short <- s$short

					m <- simple_match(sub_long, sub_short)
					m$pos_long <- c(i, i + len_short - 1)
					m$pos_short <- c(1, len_short)

					if ( self.dimer | hetero.dimer ) {
						m <- af$cal_contig_free_energy(m)
						m <- af$revert_complement(m)
					}

					return(m)
				}
				, mc.cores = af$ncores
			)

			return(lmatch)
		},
		end = function(s1, s2) {
			# Provided two strings, finds the best alignments of the shortest
			# one in the end of the longest and returns match score, sequence
			# and position
			
			# Identify long/short string
			s <- af$assign(s1, s2)
			len_short <- nchar(s$short)
			len_long <- nchar(s$long)

			# Slide substrings of the short seq on the end of the long one
			lmatch <- mclapply(1:(len_short - 1),
				FUN=function(i) {
					sub_long <- substr(s$long, len_long - i + 1, len_long)
					sub_short <- substr(s$short, 1, i)

					m <- simple_match(sub_long, sub_short)
					m$pos_long <- c(len_long - i + 1, len_long)
					m$pos_short <- c(1, i)

					if ( self.dimer | hetero.dimer ) {
						m <- af$cal_contig_free_energy(m)
						m <- af$revert_complement(m)
					}

					return(m)
				}
				, mc.cores = af$ncores
			)

			return(lmatch)
		}
	)

	# Trigger alignment
	af <- af$align()

	# Assign class
	class(af) <- 'align_fast'

	# Output
	return(af)
}

print.align_fast = function(af) print(af$matches)

# dimerization_suite -----------------------------------------------------------

# Contains functions to calculate the free energy of dimerization
dimerization_suite = function(dg_table = NULL, ncores = 1) {
	
	# Generate class
	ds <- list(
		ncores = ncores,

		get_dimers = function(s) {
			# Produce list of dimers in S
			if ( 2 > nchar(s) ) return(s)
			return(unlist(lapply(1:(nchar(s) - 1),
				FUN=function(i) substr(s, i, i + 1))))
		},

		get_dimer_free_energy = function(dimer, dgs = ds$dgs) {
			ids <- which(ds$dg_table$dimer == dimer)
			if ( 0 == length(ids) ) {
				cat(paste0('[get_dimer_free_energy] ERROR: ',
					'unavailable dG0 for current dimer ("', dimer, '")\n'))
				stop()
			}
			ds$dg_table$dG0[ids]
		},

		get_free_energy = function(seq) {
			# Returns the free energy of perfect dimerization
			# (with reverse complement)
			dgs <- unlist(lapply(ds$get_dimers(seq),
				FUN=ds$get_dimer_free_energy
			))
			return(sum(dgs))
		},

		get_self_dimer_free_energy = function(seq, only.top = F) {
			# Identifies best alignment.
			# Get longest contig.
			# Calculate free energy of the contig.
			
			dimer <- align_fast(seq, self.dimer = T, ncores = ds$ncores)
			ordered <- order.simple_match_list(dimer$matches, by='free_energy')
			dimer$matches <- dimer$matches[ordered]

			if ( only.top )
				return(min(dimer$matches[[length(dimer$matches)]]$contig$fe,
					na.rm=T))

			return(dimer)
		},

		get_hetero_dimer_free_energy = function(s1, s2, only.top = F) {
			# Identifies best alignment.
			# Get longest contig.
			# Calculate free energy of the contig.
			
			dimer <- align_fast(s1, s2, hetero.dimer=T, ncores = ds$ncores)
			ordered <- order.simple_match_list(dimer$matches, by='free_energy')
			dimer$matches <- dimer$matches[ordered]

			if ( only.top )
				return(min(dimer$matches[[length(dimer$matches)]]$contig$fe,
					na.rm=T))

			return(dimer)
		}
	)

	# Setup shortcuts
	ds$gdfe = ds$get_dimer_free_energy
	ds$gfe = ds$get_free_energy
	ds$gsdfe = ds$get_self_dimer_free_energy
	ds$ghdfe = ds$get_hetero_dimer_free_energy

	# Define standard dimerization dG table
	if ( is.null(dg_table) ) {
		ds$dg_table <- data.frame(
			dimer = c('TA', 'AT', 'AA', 'TT', 'GC', 'CG', 'GG', 'CC', 'CA',
				'AC', 'GT', 'TG', 'CT', 'TC', 'GA', 'AG'),
			dG0 = c(-0.961265, -1.474215, -1.944400, -1.944400, -3.139395,
				-3.611430, -3.069210, -3.069210, -1.953865, -1.342005,
				-1.342005, -1.953865, -1.598480, -1.574975, -1.574975,
				-1.598480)
		)
	} else {
		ds$dg_table <- read.delim(dg_table, as.is=TF)
	}

	# Assign class
	class(ds) <- "dimerization_suite"

	# Output
	return(ds)
}

################################################################################
