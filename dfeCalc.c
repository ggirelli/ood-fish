#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define ISCOMPL(A, B) ((A+B == 1) || (A+B == 5))

// Number of seq printed in the log
#define SEQ_HEAD_N 4

// Compilation: gcc -O3 dfeCalc.c -lm -o dfeCalc
// en.wikipedia.org/wiki/Nucleic_acid_thermodynamics#Nearest-neighbor_method

/**
 * Contains the dimer free energies, initialized by initFEtable.
 */
float * FEtable;

/**
 * Structure containing the settings of the current run.
 */
struct runSettings {
	// Oneline sequence file
	char * seqFile;
	// Output file
	char * outFile;
	// Sequence length
	int seqLen;
	// Energy calculation mode
	// 1: only against reverse
	// 2: only against reverse complement
	// 3: both reverse and reverse complement
	int mode;
	// Number of sequences
	int seqN;
	// Energy verbosity
	int everbose;
};

/**
 * Convert Numeric to character Base.
 * @param  b Numeric-coded DNA base.
 * @return   The original DNA base.
 */
char n2b(uint8_t b) {
	switch(b) {
		case 0:
			return 'A';
		case 1: 
			return 'T';
		case 2:
			return 'C';
		case 3:
			return 'G';
		default:
			return 'X';
	}
};

/**
 * Convert character Base to Numeric.
 * @param  b DNA base.
 * @return   Numeric-coded DNA base.
 */
int b2n(uint8_t b) {
	switch(b) {
		case 'A':
			return 0;
		case 'T':
			return 1;
		case 'C':
			return 2;
		case 'G':
			return 3;
		default:
			// Invalid character
			return 4;
	}
};

/**
 * Convert Numeric to character Sequence.
 * @param  b Numeric-coded DNA sequence.
 * @return   The original DNA sequence.
 */
uint8_t* n2s(uint8_t* seq, int seqLen) {
	for (int kk = 0; kk < seqLen; kk++)
		seq[kk] = n2b(seq[kk]);
	return seq;
}

/**
 * Convert character Sequence to Numeric.
 * @param  b DNA sequence.
 * @return   Numeric-coded DNA sequence.
 */
uint8_t* s2n(uint8_t* seq, int seqLen) {
	for (int kk = 0; kk < seqLen; kk++)
		seq[kk] = b2n(seq[kk]);
	return seq;
}

/**
 * @param  b Numeric-coded DNA base.
 * @return   Complementary base.
 */
uint8_t complementary_base(uint8_t b) {
	switch(b) {
		case 0:
			return 1;
		case 1: 
			return 0;
		case 2:
			return 3;
		case 3:
			return 2;
		default:
			return -1;
	}
};

/**
 * Log the non-0 cells of the dimer free-energy table table.
 * @param FEtable Free energy table.
 */
void printFEtable(float *  FEtable) {
	printf("\ndG table:\n");

	for (int aa = 0; aa < 4; aa++) {
		for (int bb = 0; bb < 4; bb++) {

			float dg = FEtable[ aa*1 + bb*4 ];
			if ( dg != 0 ) {
				printf(" %c%c: %f\n", n2b(aa), n2b(bb), dg);
			}

		}
	}

	printf("\n\n");
};

/**
 * Set up the free energy table (FEtable).
 */
void initFEtable() {
	// Init variable
	FEtable = (float *) malloc(16*sizeof(float));

	// Default to 0
	for (int kk = 0; kk<16; kk++)
		FEtable[kk] = 0;

	// Do it on numeric-coded bases.
	// 0 - A
	// 1 - T
	// 2 - C
	// 3 - G
	
	printf("Unit: kcal/mol\n");

	// Fill with dimer free energies
	FEtable[ 0*1 + 0*4 ] = -1.944400; // AA 
	FEtable[ 0*1 + 1*4 ] = -1.474215; // AT
	FEtable[ 0*1 + 2*4 ] = -1.342005; // AC
	FEtable[ 0*1 + 3*4 ] = -1.598480; // AG
	FEtable[ 1*1 + 0*4 ] = -0.961265; // TA
	FEtable[ 1*1 + 1*4 ] = -1.944400; // TT
	FEtable[ 1*1 + 2*4 ] = -1.574975; // TC
	FEtable[ 1*1 + 3*4 ] = -1.953865; // TG
	FEtable[ 2*1 + 0*4 ] = -1.953865; // CA 
	FEtable[ 2*1 + 1*4 ] = -1.598480; // CT
	FEtable[ 2*1 + 2*4 ] = -3.069210; // CC 
	FEtable[ 2*1 + 3*4 ] = -3.611430; // CG
	FEtable[ 3*1 + 0*4 ] = -1.574975; // GA
	FEtable[ 3*1 + 1*4 ] = -1.342005; // GT
	FEtable[ 3*1 + 2*4 ] = -3.139395; // GC
	FEtable[ 3*1 + 3*4 ] = -3.069210; // GG

	printFEtable(FEtable);
};

void printSeq(uint8_t * seq, int seqLen,  int numeric) {
	for (int kk = 0; kk < seqLen; kk++)
		if ( numeric )
			printf("%d", seq[kk]);
		else
			printf("%c", seq[kk]);
};

/**
 * Print the first SEQ_HEAD_N sequences in the input file.
 * @param seq Oneline sequence string.
 * @param numeric Whether the sequence is numeric-coded
 * @param set Run settings.
 */
void seqHead(uint8_t * seq, const struct runSettings set, int numeric) {
	printf("Sequences:\n");

	int n_seq_to_print = MIN(set.seqN, SEQ_HEAD_N);
	for (int ll = 1; ll <= n_seq_to_print; ll++) {
		printSeq(seq + ((ll-1) * set.seqLen), set.seqLen, numeric);
		printf("\n");
	}

	printf("... \n\n");
};

/**
 * Calculates the lowest dimerization free energy of two sequences
 * A and B of length seqLen. Needs also FEtable to be available.
 *
 * It navigates through the two sequences and:
 * 		- identifies complementary stretches
 * 		- calculates free energy of every stretch
 * 		- returns the free energy of dimerization of the longest stretch
 * 		  [actually, either of the longest or of the strongest]
 * 
 * @param  A 1st sequence
 * @param  B 2nd sequence
 * @param  seqLen Sequence length
 * @return   Lowest dimerization dG
 */
float energy(uint8_t * A, uint8_t *B, int seqLen, int everbose) {
	if ( everbose ) { // While DEBUGging energy calculation
		printf("\n\t");

		// Allocate string version
		uint8_t* As = (uint8_t*) malloc(seqLen*sizeof(uint8_t));
		uint8_t* Bs = (uint8_t*) malloc(seqLen*sizeof(uint8_t));
		for (int i = 0; i < seqLen; ++i) {
			As[i] = A[i];
			Bs[i] = B[i];
		}

		// Convert to string
		n2s(As, seqLen);
		n2s(Bs, seqLen);

		// Print string
		printSeq(As, seqLen, 0);
		printf("\n\t");
		for (int kk = 0; kk < seqLen; kk++) {
			if ( ISCOMPL(A[kk], B[kk]) ) {
				printf("|");
			} else {
				printf(" ");
			}
		}
		printf("\n\t");
		printSeq(Bs, seqLen, 0);
		printf("\n\n\t");
	}

	// Complementary status (whether inside a complementary stretch)
	int inC = 0;

	// Longest complementary stretch
	int best_len = 0;
	float best_emin = 0;

	// Current complementary stretch
	int cur_len = 0;
	float cur_emin = 0;

	for (int kk = 0; kk < seqLen; kk++) {

		if ( !inC ) {
			// If not in complementary stretch
			// check if entered one
			if ( ISCOMPL(A[kk], B[kk]) ) {
				inC = 1;	// Set status as inside
				cur_len++;	// Increase stretch length
			}
		} else {
			// If in a complementary stretch
			// check if it left
			inC = ISCOMPL(A[kk], B[kk]);

			if ( inC ) {
				// If still in the stretch
				// increase stretch len and dimerization free energy
				cur_len++;
				cur_emin  = cur_emin + FEtable[A[kk-1] + 4*A[kk]];
			}
			if ( !inC || kk == seqLen-1) {
				// If the compl strand ended, compare with best stretch
				if ( cur_len == best_len && cur_emin > best_emin ) {
					best_emin = cur_emin;
				}
				if ( cur_len > best_len ) {
					best_len = cur_len;
					best_emin = cur_emin;
				}

				if ( everbose ) // To DEBUG energy calculation
					printf("\n\tCur: %d - %f\n\tBest: %d - %f\n\n\t",
						cur_len, cur_emin, best_len, best_emin);

				// And leave
				cur_len = 0;
				cur_emin = 0;
			}
		}
	}

	if ( everbose ) // While DEBUGging energy calculation
		printf("\n\tFree energy:%f\n\tMax contig len:%d\n\t",
			best_emin, best_len);

	return best_emin;
};

/**
 * Calculate the max/min dimerization free energy
 * between two A and B sequences of length seqLen
 * @param  A      1st sequence
 * @param  B      2nd sequence
 * @param  seqLen Sequences length
 * @return        Free energy of dimerization
 */
float min_energy(uint8_t * A, uint8_t * B, uint64_t seqLen, int everbose) {
	// Transpose to make c1 polymers
	if ( everbose ) { // While DEBUGging energy calculation
		printf("\n---\n\n");

		// Allocate string version
		uint8_t* As = (uint8_t*) malloc(seqLen*sizeof(uint8_t));
		uint8_t* Bs = (uint8_t*) malloc(seqLen*sizeof(uint8_t));
		for (int i = 0; i < seqLen; ++i) {
			As[i] = A[i];
			Bs[i] = B[i];
		}

		// Convert to string
		n2s(As, seqLen);
		n2s(Bs, seqLen);

		// Print string
		printSeq(As, seqLen, 0);
		printf("\n");
		printSeq(Bs, seqLen, 0);
		printf("\n");
	}

	// Do the free energy calculations here
	float emin = 0;

	// A to the left of B
	if ( everbose ) printf("\n\t!!! LEFT alignment.\n");
	for (int kk = 1; kk < seqLen; kk++) {
		emin = MIN(emin, energy(A+kk, B, seqLen-kk, everbose)); 
	}

	// Centered
	if ( everbose ) printf("\n\t!!! CENTRAL alignment.\n");
	emin = MIN(emin, energy(A, B, seqLen, everbose));

	// A to the right of B
	if ( everbose ) printf("\n\t!!! RIGHT alignment.\n");
	for (int kk = 1; kk < seqLen; kk++) {
		emin = MIN(emin, energy(A, B+kk, seqLen-kk, everbose)); 
	}

	return emin;
};

/**
 * Runs the sel/hetero-dimerization calculations.
 * @param  s Run settings.
 * @return   main() output
 */
int mer(const struct runSettings s) {

	// INIT --------------------------------------------------------------------
	printf("Loading dimer free energy table...\n");
	initFEtable();

	// Sequence length
	uint64_t L = s.seqLen;

	// Oneline sequence filename
	char * fname = s.seqFile;

	// Read the sequence file to T
	printf("Reading %s\n", fname);

	// Open file pointer
	FILE * f = fopen(fname, "r");

	// Check file size
	fseek(f, 0, SEEK_END);
	uint64_t fSize = ftell(f);
	rewind(f);

	// Count sequences
	uint64_t N = fSize/L;

	// Read up to the user-defined last sequence or the last in the file.
	printf("Reading %lu sequences of length %lu...\n\n", N, L);

	// Three sequence tables:
	// 	- sequence (T)
	// 	- reverse (TR)
	// 	- complement (TC)
	// 	- reverse-complement (TRC)
	uint8_t * T = (uint8_t *) malloc(N*L*sizeof(uint8_t));
	uint8_t * TR = (uint8_t *) malloc(N*L*sizeof(uint8_t));
	uint8_t * TC = (uint8_t *) malloc(N*L*sizeof(uint8_t));
	uint8_t * TRC = (uint8_t *) malloc(N*L*sizeof(uint8_t));

	// READ AND CODE -----------------------------------------------------------

	// Read the sequence file
	size_t nread = fread(T, N*L*sizeof(uint8_t), 1, f);
	printf("Read %lu blocks of %lu bytes.\n", nread, N*L);
	fclose(f);

	// Print the first SEQ_HEAD_N sequences
	seqHead(T, s, 0);

	// Prepare the sequences for faster lookup later
	// A -> 1, T -> 2, C -> 4, G -> 8 
	// I.e., use one bit for each one, in that way the combinations of two
	// codes will be unique and can be used for hashing
	// i.e. A+A = 2, A+T = 3, T+T = 4, A+C = 5, T+C=6. ...
	printf("Translating ATCG into base2...\n");
	for (int kk = 0; kk < N*L; kk++) {
		T[kk] = b2n(T[kk]);

		// Kill process if invalid character is present.
		if ( 4 == T[kk] ) {
			printf("Invalid character read!\n");
			return (1);
		}
	}
	printf("\n");

	// Print the first SEQ_HEAD_N sequences (now binary)
	seqHead(T, s, 1);

	printf("Creating reverse strands (TR)...\n");
	for (size_t kk = 0; kk < N; kk++)
		for (int ll = 0; ll < L; ll++)
			TR[kk*L+L-ll-1] = T[kk*L+ll];

	printf("Creating complementary strands (TC)...\n");
	for (size_t kk = 0; kk < N; kk++)
		for (int ll = 0; ll < L; ll++)
			TC[kk*L+ll] = complementary_base(T[kk*L+ll]);
	
	printf("Creating reverse complementary strands (TRC)...\n");
	for (size_t kk = 0; kk < N; kk++)
		for (int ll = 0; ll < L; ll++)
			TRC[kk*L+ll] = complementary_base(TR[kk*L+ll]);

	printf("Opening %s for output...\n", s.outFile);
	FILE * outfile = fopen(s.outFile, "w");
	if ( outfile == NULL ) {
		printf("Failed to open.\n");
		return 1;
	}

	// CALC --------------------------------------------------------------------

	// Dimerization energy calculation
	uint64_t nDone = 0;			// Calculation counter
	uint64_t nTot = (N*N + N)/2;	// Total number of calculations

	printf("\nCalculating binding energies\n");

	for (uint64_t mm = 0; mm < N; mm++) { //
		for (uint64_t nn = mm; nn < N; nn++) {
			// Proper indexing would mean to take care that memory is
			// allocated only for half of the matrix (lower left or upper right)
			// A shift operator or mask has also to be used to pack the
			// single bits into each byte of E.
			
			// Default value
			float eMin = 1000;

			if ( s.mode == 1 || s.mode == 3 ) {
				// Calculate against reverse
				float ER = min_energy(T+mm*s.seqLen, TR+nn*s.seqLen,
					s.seqLen, s.everbose);
				eMin = MIN(eMin, ER);
				if ( s.everbose ) printf("\nBest: %f kcal/mol\n", eMin);
			}

			if ( s.mode == 2 || s.mode == 3 ) {
				// Calculate against complement
				float EC = min_energy(T+mm*s.seqLen, TC+nn*s.seqLen,
					s.seqLen, s.everbose);
				eMin = MIN(eMin, EC);
				if ( s.everbose ) printf("\nBest: %f kcal/mol\n", eMin);
			}

			// Print free energy to output file
			fprintf(outfile, "%lu %lu %f\n", mm+1, nn+1, eMin);
			
			// Increase counter
			nDone++;
		}

		if ( mm%1 == 0 ) {
			// Show progress
			printf("\r%3.3f %% (%lu/%lu)",
				(100.000*nDone)/nTot, nDone, nTot);
			fflush(stdout);
		}

	}
	printf("\r %3.0f %% (%lu/%lu)", floor((100.0*nDone)/nTot), nDone, nTot);
	printf("\n");

	// Close pointer to output file
	fclose(outfile);

	printf("\n ~ END ~ \n");

	free(T); //free(TR);
	free(TC); //free(TRC);
	free(FEtable);
};

int main(int argc, char ** argv) {
	// Read arguments in run settings object
	struct runSettings s;

	// Default values
	s.seqLen = 20;
	s.mode = -1;
	s.everbose = 0;

	// Read input options
	int opt = 0;
	while ((opt = getopt(argc, argv, "i:o:l:m:v::h::")) != -1) {
		switch(opt) {
			case 'h':
				printf("\nUsage: dfeCalc [-h] -i input -o output -m mode%s\n\n",
					" [-l seqLen][-v]");
				printf("Mandatory arguments:\n%s%s%s\n\n",
					"-i\tInput file\n",
					"-o\tOutput file\n",
					"-m\tMode. 1: reverse, 2: complement, 3: 1+2");
				printf("Optional arguments:\n%s%s%s\n\n",
					"-l\tSequence length, default: 20\n",
					"-v\tVerbose mode.\n",
					"-h\tShow this help page");
				return 0;
				break;
			case 'i':
				s.seqFile = optarg;
				break;
			case 'o':
				s.outFile = optarg;
				break;
			case 'm':
				s.mode = atol(optarg);
				if ( 1 > s.mode || 3 < s.mode ) {
					printf("\nInvalid -m option.\n");
					return 1;
				}
				break;
			case 'l':
				s.seqLen = atol(optarg);
				if ( 0 > s.seqLen ) {
					printf("\nInvalid -l option.\n");
					return 1;
				}
				break;
			case 'v':
				s.everbose = 1;
				break;
			case '?':
				printf("\nInvalid option received");
				break;
		}
	}

	// Check for missing options
	if ( NULL == s.seqFile ) {
		printf("Usage: dfeCalc [-h] -i input -o output -m mode %s\n",
			"[-l seqLen][-v]");
		printf("Missing mandatory input option -i.\n");
		return 1;
	}
	if ( NULL == s.outFile ) {
		printf("Usage: dfeCalc [-h] -i input -o output -m mode %s\n",
			"[-l seqLen][-v]");
		printf("Missing mandatory output option -o.\n");
		return 1;
	}
	if ( 0 > s.mode ) {
		printf("Usage: dfeCalc [-h] -i input -o output -m mode %s\n",
			"[-l seqLen][-v]");
		printf("Missing mandatory mode option -m.\n");
		return 1;
	}

	// Count sequences
	FILE * f = fopen(s.seqFile, "r");
	fseek(f, 0, SEEK_END);
	uint64_t fSize = ftell(f);
	s.seqN = fSize/s.seqLen;
	fclose(f);

	// Show current settings
	printf("\nInput file: %s", s.seqFile);
	printf("\nOutput file: %s", s.outFile);
	printf("\nSequence length: %d", s.seqLen);
	printf("\nFound %d sequences.", s.seqN);
	printf("\nMode: %d", s.mode);
	printf("\nVerbose energy calculation: %s,",
		s.everbose >= 1 ? "true" : "false");
	printf("\n\n");

	// Run the algorithm
	mer(s);
	
	printf("\n");
	return 0;
};
