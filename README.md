BCSuite
===

BCSuite contains a series of scripts, that constitute a DNA BarCode filter/design pipeline.

This document explains in detail the current pipeline for barcode filter/design. The barcodes are selected based on three features: reference genome homology, self-dimerization free energy, and hetero-dimerization free energy.

## Pipeline

1. Sequence generation.
2. [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat?command=start) against reference genome.
3. Merge BLAT output into single file.
4. Split merged BLAT output into chunks.
5. Convert chunks into blat recap output.
6. Calculate self-dimerization free energy (sdfe) for every barcode candidate.
7. Filter based on maximum homology and sdfe.
8. Calculate pair-wise hetero-dimerization free energy (hdfe) for the barcode candidates selected from the filter (7).
9. Filter based on hdfe.
10. Select barcode set with [PMC algorithm](https://www.cs.purdue.edu/homes/dgleich/codes/maxcliques/README.html).

### 1. Sequence generation

Of course, we need to produce a set of barcode candidate to be filtered and properly selected. This set can be generated in a variety of ways. One example, would be selecting all the 20 characters long substring of a set of 25-mers.

Specifically, Xu et al. published a list of 240,000 25-mers, that can be used to generate a list of 1440000 20-mers.

### 2. [BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat?command=start)

Depending on the application, it might be important for the barcodes to be genome orthogonal, i.e., being absent from the genome. In other words, the homology between a barcode and the genome should be lower than a given threshold (e.g., 60%). This can be achieved using BLAST to align each barcode to the reference genome, and then analyze the output.

Specifically, if the barcode is short (i.e., ~20 nt), the following parameters should be used:

```
blastn -query input.fa -db refDB -word_size 6 -evalue 10000 -penalty -2 -reward 1 -task 'blastn' -outfmt 4 -out blast_output.txt -num_threads X
```

### 3. Merge

To convert the output of BLAT in a more readable format (here addressed as *recap*), we first need to merge all the single barcode PSL files into a single file.

Usually, this can be easily achieved with a simple command line:

```
cat psl/* > merged_psl.psl
```

But, when the number of single barcode PSL file is too large, it is necessary to use the `psl_merge.sh`. By going in the folder containing the PSL files, just run:

```
ls | xargs -n 1000 -P 1 psl_merge.sh
```

This will generate the `../merged_psl.psl` file.

### 4. Split

Depending on the amount of available RAM on your machine, it might be difficult to operate on the `merged_psl.psl` file when many barcodes are being analyzed at the same time. We suggest then to split the file into smaller chunks using the following command line:

```
split -l 2GB  merged_psl.psl split/psl.
```

This command will split the `merged_psl.psl` file into 2GB-sized chunks of name `psl.aa`, `psl.ab`,... in the `split/` folder.

### 5. Recap

The *recap* format is useful to understand the output of the BLAT tool. Specifically, it contains one barcode (sequence) per row with the following columns:

* **i**: barcode ID number.
* **qName**: barcode ID name.
* **k**: barcode length in nt.
* **seq**: barcode sequence.
* **hX**: number of genomic loci with X homology with the barcode.
* **maxHomology**: highest homology with a genomic locus (highest X with `hX != 0`).

The script that converts the merged PSL file into the *recap* format is called `recap_blat.R`.

```
usage: recap_blat.R [--help] seq_dir psl_dir out_dir
       [--seq SEQ] [--psl PSL] [--out OUT] [--nst NST]
       [--ncores NCORES] [--suff SUFF]

Recaps BLAT output.

positional arguments:
  seq_dir           Directory containing the sequence file.
  psl_dir           Directory containing the psl file.
  out_dir           Output directory.

flags:
  -h, --help        show this help message and exit

optional arguments:
  -s, --seq SEQ     Sequence file name. [default: bc.noHeader.fa]
  -o, --out OUT     Output file name. [default: merged_table.tsv]
  -p, --psl  PSL    PSL file name. [default: merged_psl.psl]
  -n, --nst NST     Number of sub tables. [default: 1000]
  --ncores NCORES   Number of cores for parallelization. [default: 30]
  --suff SUFF       Output suffix. [default: ]
```

##### Recap after splitting

If the split step (4) was performed, the `recap_blat.R` script function should be run on every split PSL file. This can be easily achieved by running the `recap_split.sh` script.

```
usage: recap_split.sh -i psl_dir -o out_dir -s seq_dir
       [-f seq_file][-t threads][-p path]
 
Run the recap_blat.R script for every split PSL file
in the specified psl (-i) directory.

Mandatory arguments:
 -i psl_dir    Directory containing the split PSL files.
 -s seq_dir    Directory containing the seq_file (-f).
 -o out_dir    Output directory.

Optional arguments:
 -f seq_file   Name of the sequence file in seq_dir (-s).
 -t threads    Number of threads for parallelization.
 -p path Path  to recap_blat.R script. Default './recap_blat.R'
```

Then, merge the converted split recap tables with the `recap_merge.R` script.

```
usage: recap_merge.R [--help] recap_dir out_dir [--out OUT]

Merge BLAT recaps.

positional arguments:
  recap_dir         Directory containing the recap files.
  out_dir           Output directory.

flags:
  -h, --help        show this help message and exit

optional arguments:
  -o, --out OUT     Output file name. [default: merged_table.tsv]
```

### 6. Calculate sdfe

If needed, it is also possible to add a column with the calculated self-dimerization free energy (sdfe) to the recap table, using the `sdfe_calc.R` script.

```
usage: sdfe_calc.R [--help] rec_dir out_dir
       [--rec REC] [--out OUT] [--ncores NCORES] [--suff SUFF]

Calculate self dimerization free energy.

positional arguments:
  rec_dir               Directory containing the recap table.
  out_dir               Output directory.

flags:
  -h, --help            show this help message and exit

optional arguments:
  -r, --rec REC         Recap table file name. [default: merged_table.tsv]
  -o, --out OUT         Output file name. [default: merged_table.sdfe.tsv]
  -n, --ncores NCORES   Number of cores for parallelization. [default: 30]
  -s, --suff SUFF       Output suffix. [default: ]
```

### 7. Filter (maxHomology & sdfe)

Then, the recap table with sdfe column should be filtered based on sdfe and maxHomology with the `sdfe_filter.R` script. Use the output of the previous step (6) as the input (fin). The default SDFE threshold is set at -5 kcal/mol while the maxHomology at 75% of the sequence.

```
usage: sdfe_filter.R [--help] in_dir fin out_dir
       [--out OUT] [--homo_thr HOMO_THR] [--sdfe_thr SDFE_THR] [--ncores NCORES] [--seqlib SEQLIB]

Filter BLAT recap table based on SDFE and maxHomology.

positional arguments:
  in_dir          Input directory.
  out_dir         Output directory.

flags:
  -h, --help      show this help message and exit

optional arguments:
  -f, --fin FIN   Input file name.
                  [default: merged_table.sdfe.tsv]
  -o, --out OUT   Output file name.
                  [default: merged_table.filtered.sdfe.tsv]
  --ht HT         Homology threshold.
                  [default: 0.75]
  -s, --st ST     Self-dimerization free energy threshold (positive).
                  [default: -5]
  --sl SL         Path to seq_lib.R
                  [default: ./seq_lib.R]
```

### 8. Calculate hdfe

To calculate the hetero-dimerization free energy (hdfe) of every pair of filtered barcode candidate, use the `dfeCalc` script (based on code by Dr. Erik Wernersson).

```
usage: dfeCalc [-h] -i input -o output -m mode [-l seqLen][-v]

Mandatory arguments:
-i  Input file
-o  Output file
-m  Mode. 1: reverse, 2: complement, 3: 1+2

Optional arguments:
-l  Sequence length, default: 20
-v  Verbose mode.
-h  Show this help page

```

Using the `-v` option to trigger verbose mode and have more details on the free energy calculations. The three different modes are:

1. Calculate the HDFE of every oligo and the reverse of every other.
2. Calculate the HDFE of every oligo and the reverse complement of every other.
3. Keep the minimum (strongest bind) HDFE from the other two modes.

The script **input file** is a oneline fasta file with the sequences from the filtered table generated at step 7. This file can be produced using the following command line:

```
cat merged_table.filtered.sdfe.tsv | cut -f 4 | sed 1d | tr -d "\n" > sequence.filtered.oneline.fa
```

### 9. Filter (hdfe)

Filter the HDFE table using the `hdfe_filter.R` script. The default threshold is set at -9 kcal/mol.

```
usage: hdfe_filter.R [--help] in_dir out_dir
       [--fin FIN][--out OUT][--pmc PMC][--ht HT]

Filter BLAT recap table based on SDFE and maxHomology.

positional arguments:
  in_dir            Input directory.
  out_dir           Output directory.

flags:
  -h, --help        show this help message and exit

optional arguments:
  -f, --fin FIN     Input file name.
                    [default: hdfe.dat]
  -o, --out OUT     Output file name.
                    [default: hdfe_filtered.dat]
  -p, --pmc PMC     Output PMC file name.
                    [default: hdfe_filtered_graph.dat]
  --ht HT           Homology threshold (positive).
                    [default: -9]
```

### 10. Select barcode

Finally, select a set of compatible barcodes with the [`pmc` algorithm](https://www.cs.purdue.edu/homes/dgleich/codes/maxcliques/README.html).

```
./pmc -f hdfe_filtered_graph.dat -t ncores > pmc.log
```

```
Usage: ./pmc -a alg -f graphfile -t threads -o ordering -h heu_strat -u upper_bound -l lower_bound -r reduce_wait_time -w time_limit 
  -a algorithm                 : Algorithm for solving MAX-CLIQUE: 0 = full, 1 = no neighborhood cores, 2 = only basic k-core pruning steps  
  -f graph file                : Input GRAPH file for computing the Maximum Clique (matrix market format or simple edge list). 
  -o vertex search ordering    : Order in which vertices are searched (default = deg, [kcore, dual_deg, dual_kcore, kcore_deg, rand]) 
  -d decreasing order          : Search vertices in DECREASING order. Note if '-d' is not set, then vertices are searched in increasing order by default. 
  -e neigh/edge ordering       : Ordering of neighbors/edges (default = deg, [kcore, dual_deg, dual_kcore, kcore_deg, rand]) 
  -h heuristic strategy        : Strategy for HEURISTIC method (default = kcore, [deg, dual_deg, dual_kcore, rand, 0 = skip heuristic]) 
  -u upper_bound               : UPPER-BOUND on clique size (default = K-cores).
  -l lower_bound               : LOWER-BOUND on clique size (default = Estimate using the Fast Heuristic). 
  -t threads                   : Number of THREADS for the algorithm to use (default = 1). 
  -r reduce_wait               : Number of SECONDS to wait before inducing the graph based on the unpruned vertices (default = 4 seconds). 
  -w time_limit                : Execution TIME LIMIT spent searching for max clique (default = 7 days) 
  -k clique size               : Solve K-CLIQUE problem: find clique of size k if it exists. Parameterized to be fast. 
  -s stats                     : Compute BOUNDS and other fast graph stats 
  -v verbose                   : Output additional details to the screen. 
  -? options                   : Print out this help menu.
```
