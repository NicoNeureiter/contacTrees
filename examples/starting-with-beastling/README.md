# Adding contacTrees to a BEAST xml

In this example, I create a minimal XML file for a linguistic phylogeny using BEASTling on readily available data. The XML is constructed with per-concept rate variation and otherwise minimal, to ensure that the XML file contains separate alignments for different concepts (which is necessary for adding contacTrees) but is otherwise relatively easy to work with.

## Preparation

The best way to use contacTrees is to add the necessary BEAST2 objects to an existing configuration file. The focus of this example is on the contacTrees, not on the other elements of configuration, so I use a standard dataset (Grollemund et al.'s Bantu data) with a standard tool (BEASTling) to generate a base XML, which I will then extend to be a contacTrees model below.

### Data

Of the four large language families with extensive phylogenetic studies, Bantu stands out for haviing an easily available data set (Indo-European's picture is much more complicated) which does not include forms from other language families (like the big Austronesian dataset, ABVD, does) and dialect chains seem to be a lesser factor in its formation (as opposed to Sino-Tibetan).

I therefore base this exapmle on the Bantu dataset from

> Grollemund, Rebecca & Branford, Simon & Bostoen, Koen & Meade, Andrew & Venditti, Chris & Pagel, Mark. 2015. Bantu expansion shows that habitat alters the route and pace of human dispersals. Proceedings of the National Academy of Sciences 112(43). 13296–13301. (doi:10.1073/pnas.1503793112)

in the CLDF version available from https://github.com/lexibank/grollemundbantu ([v1.0rc6](https://github.com/lexibank/grollemundbantu/releases/tag/v1.0rc6)). The data set is available under the CC-BY-NC 4.0 International license, so I have included a copy of it in [this folder](cldf/).

### BEASTling

To create the XML file, I use the python helper [BEASTling](https://github.com/lmaurits/BEASTling/), which from a small description file can generate a full XML file when given a CLDF dataset.

> Maurits, Luke & Forkel, Robert & Kaiping, Gereon A. & Atkinson, Quentin D. 2017. BEASTling: a software tool for linguistic phylogenetics using BEAST 2. PLOS ONE 12(8). e0180908. (doi:10.1371/journal.pone.0180908)

I use the most recent published version of BEASTling, v1.5.1 (installed using `pip install beastling`) with the [following configuration file](bantu-beastling.conf).

> ```
> [admin]
> basename=bantu-beastling
> [model bantu-vocabulary]
> model = binaryctmc
> data = cldf/cldf-metadata.json
> rate_variation=True
> ```
> — bantu-beastling.conf

Generating the BEAST2 XML file using `beastling bantu-beastling.conf` creates [the basic XML file in this folder](bantu-beastling.xml), which will serve as the basis for adding contacTrees in the following steps.

## Adding contacTrees

...
