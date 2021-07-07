---
title: "scCAN package manual"
author: "Bang Tran"
date: "2021-07-07"
output:
  html_document:
    df_print: paged

editor_options:
  chunk_output_type: inline
vignette: >
  %\VignetteIndexEntry{scCAN package manual}
  %\VignetteEngine{knitr::knitr}
  %\usepackage[UTF-8]{inputenc}

---

# Introduction
Clustering is the key step to define cell types from the heterogeneous cells population. One critical challenge is that single-cell RNA sequencing (scRNA-Seq) experiments generate a massive amount of data with an excessive noise level, which hinders many computational algorithms. We propose a single cell Clustering method using Autoencoder and Network fusion (scCAN) that can accurately segregate cells from the high-dimensional noisy scRNA-Seq data.

# Installation
### Install scDHA package and other requirements
scDHA can be installed fron GitHub using below instruction. scDHA depends on the `keras` package in python to build and train the autoencoders. You need to have your python environment set up before installing `keras`. You can install `miniconda` to quickly set it up using below code.


# Analysis on the sample dataset


```r
library(scCAN)
#Load example data (SCE dataset)
data("SCE")

#Get data matrix and label
data <- t(SCE$data); label <- as.character(SCE$cell_type1)
```

### Clustering

```r
#Generate clustering result, the input matrix has rows as samples and columns as genes
result <- scCAN(data, r.seed = 1)

#The clustering result can be found here 
cluster <- result$cluster

#Calculate adjusted Rand Index using mclust package
ari <- round(scCAN::adjustedRandIndex(cluster,label), 2)
print(paste0("ARI = ", ari))
```

