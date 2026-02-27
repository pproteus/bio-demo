This is a work-in-progress.

Motivation: What phylogenetic reconstructions are replicable with only a minimal amount of domain knowledge? We're not trying to make new discoveries here, rather we are trying to build lay-confidence in robust results.

In this demo: 
- download mtDNA sequences from NCBI
- extract genes from those sequences
- align those genes across species
- make a phylogenetic tree out of the alignment
- ...

Requirements:
- [biopython](https://biopython.org)
- [muscle](https://github.com/rcedgar/muscle)
- [iqtree](https://iqtree.github.io)
- [a way to visualize the treefile you just made](https://github.com/husonlab/dendroscope3)