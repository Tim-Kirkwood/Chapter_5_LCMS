# Annotation 
Annoattion was achieved using IPApy2 via **run_IPA.py**.  This generated (i) a clustered peak table where peaks were grouped according to intensity correlation, and (ii) an annotation table containing puative annotations for the peaks in the peak table.  These were then parsed for BGC-specific annotations via **search_known_annotations.py**.  Any BGC-specific annotations would be stored in ```specific_subgraphs```, but none were identified. 

Note, a seperate subfolder **"DB"** should be made in the same directory containing the **run_IPA.py** and **search_known_annotations.py**.  The DB subfolder should contain the files **"IPA_MS1.csv"**, **"adducts.csv"**, and **"allBio_connections.csv"**.  These can be downloaded from the IPApy2 repository (https://github.com/francescodc87/ipaPy2/tree/main) and are described in the paper below for IPApy2:

Del Carratore F, Eagles W, Borka J, Breitling R. ipaPy2: Integrated Probabilistic Annotation (IPA) 2.0-an improved Bayesian-based method for the annotation of LC-MS/MS untargeted metabolomics data. Bioinformatics. 2023 Jul 1;39(7):btad455. doi: 10.1093/bioinformatics/btad455. PMID: 37490466; PMCID: PMC10382385.
