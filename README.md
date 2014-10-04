biotcm-cspn
===========

A method to build pathway interaction network.

## Environment requirements

1. Linux/Unix OS, Windows OS
2. Ruby
3. R

## Inputs

1. a list of interested genes

2. PPI network [Symbol	Symbol]

3. concerned pathway list [KEGG pathway ids]

## Steps



## Outputs

1. Pathway Interaction Network
	e.g. ./result/pathway_interaction_network.txt
	format: pathway1_name	pathway2_name	p_value

2. Pathway Pairs involved in each related active PPI
	e.g. ./result/aPPIs_2_pathway_pairs.txt
	format: gene1_Entrez_ID gene2_Entrez_ID pathway_pairs_list

## FAQs

Q: How to use a specific PPI network?

A: 

Q: How to provide new concerned pathway list or GO term list?

A: 

Q: How to provide new files of gene sets corresponding to pathways or GO terms?

A:

## Reference

Yezhou Huang and Shao Li, “Detection of characteristic sub pathway 
network for angio-genesis based on the comprehensive pathway network”, 
BMC Bioinformatics, Volume 11 Supplement 1, 2010: Selected articles 
from the Eighth Asia-Pacific Bioinformatics Conference, Bangalore, 
India, January 2010 (APBC 2010).
