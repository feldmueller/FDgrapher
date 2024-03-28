# Frame and Discourse Grapher (FDgrapher)

The Frame and Discourse Grapher (FDgrapher) is a result of my PhD thesis ("Das Framing von Extremismusvarianten im medialen Diskurs der Jahre 1999 – 2021. Eine corpus-driven Methode zur Erschließung und Visualisierung semantischer Frames", handed in on 28.03.2024 at Leipzig University).

## Installation

* create a new Python environment and clone this repository
* install requirements.txt for Python (`pip install requirements.txt`)
* install the [Corpus Workbench](https://cwb.sourceforge.io) (CWB). You can find an excellent (German) tutorial [here](https://www.bubenhofer.com/korpuslinguistik/kurs/index.php?id=cwb_start.html).
* install the libraries of fdgrapher.R (`install.packages("PACKAGE")`)

## Prerequisites

* trained Word2Vec Word Embedding Models
    * in the provided example, you can train models using [TWEC](https://github.com/valedica/twec) 
* a CWB corpus
    * in the provided example, you can download the (German) GermaParl corpus [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10421773.svg)](https://doi.org/10.5281/zenodo.10421773) via polmineR [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4042093.svg)](https://doi.org/10.5281/zenodo.4042093).

## Usage

* The FDgrapher is based on Python and R.
* fdgrapher_example.ipynb and fdgrapher_example.R illustrate the use of FDgrapher, using the GermaParl corpus. You can adapt the scripts to ur needs.
    * especially, you will have to specify paths to your trained Word2Vec Word Embedding Models.

FDgrapher can be used in 5 steps with the provided example files:

1. Export the .vrt file of your CWB corpus using `cwb-decode -Cx -r PATH/registry GERMAPARL -ALL > germaparl.vrt` on the command line
2. (if you haven't trained your Word Embedding Models yet, do so now)
3. Enter parameters in the file and run `fdgrapher_example.ipynb`. Your .vrt will be enriched with annotations of word embedding clusters.
4. Reimport the .vrt file as a new CWB corpus (make sure to enter all relevant attributes) using
 

`cwb-encode -s -c utf8 -x -v -d PATH/CORPUS/germa_clusters -f path/germaparl_decades_clusters.vrt -R PATH/REGISTRY/registry/germa_clusters -P pcluster -P dcluster -V corpus -V party -V parliamentary_group  -V speaker -V lp -V session -V date -V role -V interjection -V agenda_item -V agenda_item_type -V src -V url -V decade -V year`

and then

`cwb-makeall -r PATH/cwb-3.4.22/registry GERMA_CLUSTERS`

5. Enter parameters in the file and run `fdgrapher_example.R`

## Interpretation

* Word Embedding clusters can be read as elements of Semantic Frames and – combined with a collocation analysis – model complex frame / discourse structures
* Nodes represent frame elements and edges relations between the frame elements
* Tooltips provide keyness and frequency information of the contained words
