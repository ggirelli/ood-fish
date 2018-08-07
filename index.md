Orthogonal Oligo Design for Fluorescence *In Situ* Hybridization (**OOD-FISH**) is a **pipeline** for the <b>identification of oligonucleotide sequences that do <u>not</u> hybridize to a reference genome of interest</b>. Such orthogonal oligonucleotides can easily be used as barcodes for PCR or similar protocols.

OOD-FISH can be used to identify the largest set of oligonucleotides (from the input) that:

* do not hybridize to the genome of interest,
* do not self-dimerize, and
* do not hetero-dimerize.

Please, note that OOD-FISH does not include any filters based on secondary-structures of the oligonucleotide. This can be easily achieved with available software, like [OligoArrayAux](http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux).

The pipeline includes the following steps:

1. Sequences (input) generation.
2. Alignment against reference genome.
3. Calculate self-dimerization free energy (SDFE) of every barcode candidate.
4. Filter candidates based on maximum homology and SDFE.
5. Calculate pair-wise hetero-dimerization free energy (HDFE) of selected candidates.
6. Filter based on HDFE.
7. Select largest barcode set with [Parallel Maximum Clique (PMC) algorithm](https://github.com/ryanrossi/pmc).

## Useful links

* [Install](https://ggirelli.github.io/ood-fish/install)
* [Usage](https://ggirelli.github.io/ood-fish/usage)
* [Contributing Guidelines](https://ggirelli.github.io/ood-fish/contributing)
* [Code of Conduct](https://ggirelli.github.io/ood-fish/code_of_conduct)
