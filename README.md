OOD-FISH v0.0.2
===

Orthogonal Oligo Design for Fluorescence *In Situ* Hybridization (**OOD-FISH**) is a pipeline for the identification of oligonucleotide sequences that do not hybridize to a reference genome of interest. Such orthogonal oligonucleotides can easily be used as barcodes for PCR or similar protocols.

Read the [documentation](https://ggirelli.github.io/ood-fish/) for more details.

Installation
---

1. Clone the git repository locally.

```bash
git clone https://github.com/ggirelli/ood-fish/
cd ood-fish
```

2. Install R dependencies.

```R
for ( p in c("argparser", "data.table", "parallel", "readr") )
    if ( !require(p, character.only = T) )
        install.packages(p)
```

3. Compile if needed.

```bash
gcc src/dfeCalc.c -o src/dfeCalc
```

4. Install PMC algorithm following the instructions [here](https://github.com/ryanrossi/pmc#setup).

Usage
---

More details on how to run **OOD-FISH** are available in the [documentation](https://ggirelli.github.io/ood-fish/).

Contributing
---

We welcome any contributions to `ood-fish`. Please, refer to the [contribution guidelines](https://ggirelli.github.io/ood-fish/contributing) if this is your first time contributing! Also, check out our [code of conduct](https://ggirelli.github.io/ood-fish/code_of_conduct).

License
---

```
MIT License
Copyright (c) 2017-2018 Gabriele Girelli
```
