# bundle
Bundle for viral transmission network

## Preprocessing
1) Download samples from NCBI
2a) cat them together in one txt file
2b) rename samples 
3) MSA them using MAFFT
4) Run RAXML to generate tree using MSA

## Forming of clusters: 
5) Run 1distanceMatrixFormation.py to form dist matrix from MSA 
6) Run 2MDS+cluster.py

## Evaluation (ARI values, phylotree):
7) Run 3randindex.py (please run correctPhylo() first before randindexing phylo)
8) Run 4treecolouring.R (using RAXML tree generated earlier)

## seqTrack network formation: (this is time consuming)
0) Obtain LoFreq data
1) Run 5filepartition.py to partition samples into diff folders w accession list  
2) Run 6faster.R, which gives a Dmat HellingerMatrix.Rda
3) Using the Dmat, run 7seqTrack.R to plot seqTrack network
