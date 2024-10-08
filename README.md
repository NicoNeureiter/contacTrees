# contacTrees
A [BEAST 2](http://beast2.org/) package for jointly inferring phylogies and contact events from linguistic data.

## Installation
* Install BEAST 2 (available from [http://beast2.org](http://beast2.org)).
* Install the contacTrees package through the [package manager](http://www.beast2.org/managing-packages/)

## Using contacTrees
To use contacTrees, create a BEAST2 file that has separate alignments for different concepts (eg. using [BEASTling](https://github.com/lmaurits/BEASTling/) with rate variation or by filling in an existing template using [Lexedata](https://github.com/Anaphory/lexedata)) and add operators, state nodes, and change likelihood elements. A detailed example can be found inside [the examples folder](examples/starting-with-beastling).

## Case study
Scripts and data for reproducing a case study using contacTrees is avaialbe at [github.com/NicoNeureiter/contacTrees-IndoEuropean](https://github.com/NicoNeureiter/contacTrees-IndoEuropean)

## Reference
Neureiter N, Ranacher P, Efrat-Kowalsky N, Kaiping GA, Weibel R, Widmer P, Bouckaert RR. Detecting contact in language trees: a Bayesian phylogenetic model with horizontal transfer. Humanities and Social Sciences Communications. 2022 Jun 17;9(1):1-4.[DOI:10.1057/s41599-022-01211-7](https://doi.org/10.1057/s41599-022-01211-7)
PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-1262191/v1]
