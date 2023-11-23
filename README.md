# MRggi
[![DOI](https://zenodo.org/badge/717145137.svg)](https://zenodo.org/doi/10.5281/zenodo.10108230)

MR-GGI can accurately infer gene-gene interactions using Mendelian randomization.

# Library
```
#remotes::install_github("hiows/MRggi")  # install MRggi
library(MRggi)
library(dplyr)
library(data.table)
```

# Example
1. data load
```
y.chr1 = fread("chr1_Gene.txt") %>% as.matrix()  # SNP genotype
X.chr1 = fread("chr1_SNP.txt") %>% as.matrix()  # Gene expression
```
2. Fine mapping
```
X.chr1_FineMap = FineMapping(X = X.chr1, y = y.chr1)
```
3. MR ggi
```
MRggi.results = MRggi(y = y.chr1, X = X.chr1_FineMap, cor.thr = 0.8)
```
3. Convert MRggi result for GRN
```
netInfo = convertnet(res = MRggi.results, fdr.thr = 0.05, c.size = 2)
```
4. GRN visualization
```
visGRN(nodes = netInfo$nodes, edges = netInfo$edges)
```
#  Further analysis  
1. GO enrichment analysis
2. KEGG pathway analysis
