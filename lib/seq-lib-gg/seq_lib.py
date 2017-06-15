'''
@description: Base library for sequencing data analysis.
@author: Gabriele Girelli
@version: 0.1.0
'''

# FUNCTIONS ====================================================================

def align_fast(sa, sb):
	'''
	Find best alignments of shortest sequence in the longest.

	Args:
		sa (string): first string.
		sb (string): second string.

	Return:
		list: best scrore and alignment.
	'''

	if (len(sb) > len(sa)):
		s1 = sb
		s2 = sa
	else:
		s1 = sa
		s2 = sb

	topscore = 0
	seqs = []
	for i in range(len(s1)-len(s2)+1):
		sub = s1[i:(i+len(s2))]
		score = 0.0
		for j in range(len(sub)):
			if sub[j] == s2[j]:
				score += 1
		score = score / len(s2)
		if score > topscore:
			topscore = score;
			seqs = [(sub, s2, i)]
		elif score == topscore:
			seqs.append((sub, s2, i))

	return([topscore, seqs])

def print_align(sa, sb):
	'''
	Print a readable alignment of two strings.

	Args:
		sa (string): first string.
		sb (string): second string.

	Return:
		NULL: print the alignment in a readable way.
	'''

	if (len(sa) != len(sb)):
		print 'ERROR'
		return
	else:
		c = 0
		s = sa + '\n'
		for i in range(len(sa)):
			if (sa[i] == sb[i]):
				s += '|'
				c += 1
			else:
				s += ' '
		s += '\n' + sb + '\n'
		s = 'Score: ' + str(c) + '\n' + s
		print(s)
		return s

def rc(na, t):
	'''
	Args:
		na (string): nucleic acid sequence.
		t (string): nucleic acid type, either 'dna' or 'rna'.

	Return:
		string: reverse complement of na.
	'''

	# Identify type
	t=t.lower()

	# Select alphabet
	if t == 'dna':
		ab = ["ATCG", "TAGC"]
	elif t == 'rna':
		ab = ["AUCG", "UAGC"]
	else:
		print('ERROR: unknown na type.')
		return

	rab = ab[1].strip().lower()
	ab = ab[0].strip().lower()

	# Check provided string
	na = na.lower()
	for c in na:
		if not c in ab:
			print 'ERROR: provided string conflicts with the selected alphabet.'
			return

	# Calculate reverse
	r = na[::-1]

	# Calculate reverse complement
	rc = []
	for c in r:
		rc.append(rab[ab.index(c)])
	rc=''.join([str(c) for c in rc]).upper()

	return(rc)

################################################################################
