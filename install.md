---
title: "OOD-FISH Installation"
---

# How to install OOD-FISH

* Clone the git repository locally.

```bash
git clone https://github.com/ggirelli/ood-fish/
cd ood-fish
```

* Install R dependencies by running the following in an R command line.

```R
for ( p in c("argparser", "data.table", "parallel", "readr") )
    if ( !require(p, character.only = T) )
        install.packages(p)
```

* Compile if needed.

```bash
gcc src/dfeCalc.c -o src/dfeCalc
```

* Install PMC algorithm following the instructions [here](https://github.com/ryanrossi/pmc#setup).

