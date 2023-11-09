# GOTHiC: Gene Ontology Topology from Hi-C

GOTHiC is a python tool that allows for the creation, annotation, and analysis of Hi-C networks. In essence, it leverages graph theory techniques to allow for the assessment of structure-function relationships in nuclear chromatin.

![A summary of GOTHiC's primary workflow](https://github.com/LavalleeAdamLab/GOTHiC/assets/61287366/7d13adaa-5a70-4c1b-b0a7-7116a9736b1b)

## Using GOTHiC:
To use GOTHiC, copy `gothic.py` into a directory within your python PATH, and simply 'from import gothic *'. Individual functions are available for specific tasks involving the conversion of .sam or .fastq files to Hi-C graphs, graph functional annotation with Gene Ontology terms or REACTOME metabolic pathways, clustering analysis, and visualization. The most conventional use case for GOTHiC can be invoked by calling the function gothic_full_run() with appropriate arguments for your use case.





