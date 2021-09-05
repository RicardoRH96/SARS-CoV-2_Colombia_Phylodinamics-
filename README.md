# SARS-CoV-2_Colombia_Phylodinamics-
Here you will find code, XMLs, Logs, Trees and other files generated for the Phylodinamic, molecular and phylogeographic tracking/analysis of SARS-CoV-2 pandemic in Colombia.

## Pipeline for basic  Birth Death Skyline Serial and Biceps analysis 

Run with BEAST using the -DF option with our updated alignment defined in an JSON format which has 2 demes definitions (trait2demes) [here](/templates/).


1) Run with BEAST 2.6.3pre

```
/path/to/beast/bin/beast -DF ../sequences/dataset.json ..SARS-CoV-2_Colombia_Phylodinamics-
/BDSKL2.xml
```
