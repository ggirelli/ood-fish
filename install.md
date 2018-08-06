---
title: "OOD-FISH Installation"
---

# How to install OOD-FISH

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

