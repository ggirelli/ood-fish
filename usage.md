---
title: "OOD-FISH Usage"
---

# How to run OOD-FISH

## 1. Sequence generation

OOD-FISH needs as input a set of oligonucleotides (barcode candidates). This set can be generated in a variety of ways.

One example, would be selecting all the 20 characters long substring of a set of 25-mers. Specifically, [Xu et al.](http://www.pnas.org/content/106/7/2289) published a list of 240,000 25-mers, that can be used to generate a list of 1440000 20-mers.

## 2. Align against reference genome of interest

Depending on the application, it might be important for the barcodes to be genome orthogonal, i.e., being absent from the genome. In other words, the homology between a barcode and the genome should be lower than a given threshold (e.g., 60%).

Information on the homology of each candidate against the reference genome of interest can be easily obtained running short-read aligners like [BLAT](http://genome.ucsc.edu/goldenPath/help/blatSpec.html) or [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download). OOD-FISH supports both tools.

### 2.a Using BLAST

```bash
blastn -query input.fa -db refDB -word_size 6 -evalue 10000 -penalty -2 -reward 1 -task 'blastn' -outfmt 6 -out blast_output.txt -num_threads X
```

Then, we convert the output of **BLAST** in a more readable format (here addressed as *recap*). The *recap* format is useful to understand the output of the BLAT tool. Specifically, it contains one barcode (sequence) per row with the following columns:

* **i**: barcode ID number.
* **qName**: barcode ID name.
* **k**: barcode length in nt.
* **seq**: barcode sequence.
* **hX**: number of genomic loci with X homology with the barcode.
* **maxHomology**: highest homology with a genomic locus (highest X with `hX != 0`).

The script that converts the output of BLAST (when using `-outfmt 6`) into the *recap* format is called `recap_blast.R`. Run `recap_blast.R -h` for more details on how to run the script.

### 2.b Using BLAT

```bash
blat -tileSize=6 -stepSize=1
-minMatch=1 -oneOff=1 -minScore=0 -minIdentity=0 -maxGap=0 -repMatch=131071
-noHead genome.fa input.fa blat_output.psl
```

Then, we convert the output of **BLAT** into the *recap* format. The script that converts the merged PSL file into the *recap* format is called `recap_blat.R`. Run `recap_blat.R -h` for more details on how to run the script.

#### BLAT caveats

Differently from BLAST, BLAT is not parallelized by default. Parallelization can still be easily achieved by exploiting software like [GNU Parallel](https://www.gnu.org/software/parallel/) and feeding it one oligomer at a time (batch of two lines from the input fasta).

In this case, to convert the output of BLAT into the *recap* format, we first need to merge all the single barcode PSL files into a single file.

Usually, this can be easily achieved with a simple command line:

```bash
cat psl/* > merged_psl.psl
```

But, when the number of single barcode PSL file is too large, it is necessary to use the `psl_merge.sh`. First go in the folder containing the PSL files, then run:

```bash
ls | xargs -n 1000 -P 1 psl_merge.sh
```

This will generate the `../merged_psl.psl` file.

Still, depending on the amount of available RAM on your machine, it might be difficult to operate on the `merged_psl.psl` file when many barcodes are being analyzed at the same time. We suggest then to split the file into smaller chunks using the following command line:

```bash
split -l 2GB  merged_psl.psl split/psl.
```

This command will split the `merged_psl.psl` file into 2GB-sized chunks of name `psl.aa`, `psl.ab`,... in the `split/` folder.

If splitting was necessary, the `recap_blat.R` script function should be run on every split PSL file. This can be easily achieved by running the `recap_split.sh` script.

Then, merge the converted split recap tables with the `recap_merge.R` script.

## 3. Calculate Self-Dimerization Free Energy (SDFE)

If needed, it is also possible to add a column with the calculated self-dimerization free energy (SDFE) to the recap table, using the `sdfe_calc.R` script. `sdfe_calc.R -h` for more details.

## 4. Filter based on homology and SDFE

Then, the *recap* table with SDFE column should be filtered based on SDFE and maxHomology with the `sdfe_filter.R` script. Use the output of the previous step (6) as the input (fin). The default SDFE threshold is set at -5 kcal/mol while the maxHomology at 75% of the sequence. `sdfe_filter.R -h` for more details.

## 5. Calculate pair-wise Hetero-Dimerization Free Energy (HDFE)

To calculate the hetero-dimerization free energy (HDFE) of every pair of filtered barcode candidate, use the `dfeCalc` script (based on code by Dr. Erik Wernersson). `dfeCalc -h` for more details.

Using the `-v` option triggers verbose mode and provides more details on the free energy calculations.

Three (3) different free-energy calculation modes are available:

1. Calculate the HDFE of every oligo and the reverse of every other.
2. Calculate the HDFE of every oligo and the reverse complement of every other.
3. Keep the minimum (strongest bind) HDFE from the other two modes.

The script **input file** is a oneline fasta file with the sequences from the filtered table generated at step 4. This file can be produced using the following command line:

```bash
cat table.filtered.sdfe.tsv | cut -f 4 | sed 1d | tr -d "\n" > sequence.filtered.oneline.fa
```

## 6. Filter based on HDFE

Filter the HDFE table using the `hdfe_filter.R` script. The default threshold is set at -9 kcal/mol. `hdfe_filter.R -h` for more details.

## 7. Select largest set of compatible barcodes

Finally, select a set of compatible barcodes with the [`PMC` algorithm](https://github.com/ryanrossi/pmc).

```
./pmc -f hdfe_filtered_graph.dat -t ncores > pmc.log
```
