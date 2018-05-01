# hmR

`hmR` is an R package to analyze [`heteromotility`](https://github.com/cellgeometry/heteromotility) data.  
The semantics of `hmR` are inspired by [`Seurat`](https://github.com/satijalab/seurat).

## Installation 

```
require(devtools)
install_github('jacobkimmel/hmR')
```

## Usage

`hmR` semantics focus on a `heteromotility` data set object, initialized from raw cell behavior feature data extracted using [`heteromotility`](https://github.com/cellgeometry/heteromotility).

```
library(hmR)

df = read.csv('path/to/motility_statistics.csv')
mot = hmMakeObject(raw.data=df)

# Perform hierarchical clustering
mot = hmHClust(mot, k = 3, method='ward.D2')

# Run and plot PCA
mot = hmPCA(mot)
mot = hmPlotPCA(mot)

# Run and plot tSNE
mot = hmTSNE(mot)
mot = hmPlotTSNE

# etc.
```
