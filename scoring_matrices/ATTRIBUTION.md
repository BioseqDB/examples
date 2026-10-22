# Data and reference matrices

The 1,000 sampled mitochondrial records come from Li et al., *Variation and association to diabetes in 2000 full mtDNA sequences mined from an exome study in a Danish population*, European Journal of Human Genetics 22, 1040–1045 (2014), DOI [10.1038/ejhg.2013.282](https://doi.org/10.1038/ejhg.2013.282). The manifest selects the first 1,000 records explicitly titled complete genome from the study's accession interval KF161060–KF163059. Partial-genome records passed during that selection are listed in `selection_exclusions.csv`. The tutorial retains public submitted isolate identifiers, accession versions, organism names and lengths in `samples.csv`. Sequence bodies are downloaded at runtime.

The rCRS sequence and coding annotations come from [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1). Translation uses [NCBI genetic code 2](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG2), including its start-codon convention and reverse-complement CDS on the opposite strand.

BLOSUM62 and PAM250 joint probabilities and integer scores come from the NCBI C++ Toolkit, credited in `reference_matrices/matrix-provenance.json`. The original NCBI public-domain notice is preserved at `reference_matrices/NCBI-PUBLIC-DOMAIN-NOTICE.txt`.

The BLOSUM construction was introduced by [Henikoff and Henikoff (1992)](https://doi.org/10.1073/pnas.89.22.10915). The original PAM model was introduced by Dayhoff and colleagues. The tutorial's empirical mitochondrial matrices are recalculated from the downloaded cohort using the SQL in this directory.
