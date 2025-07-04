# micro_phyloplace
Phylogenetic placement of amplicons/reads using WISE2/MAFFT/EPA_ng (Docker/micromamba stack)

This software includes data from the UniProt Knowledgebase (UniProtKB/Swiss-Prot), licensed under CC BY 4.0. © UniProt Consortium. See https://www.uniprot.org for details.
Additionally, MIBiG data is included as part of the default reference phylogeny/alignment  (citation: Terlouw et al., MIBiG 3.0: a community-driven effort to annotate experimentally validated biosynthetic gene clusters, Nucleic Acids Research, Volume 51, Issue D1, 6 January 2023, Pages D603–D610, https://doi.org/10.1093/nar/gkac1049).

# Installation

Have docker installed on your Linux system. Clone the repository. From its main directory, run the build script with:
```
docker build ./ -t micro_phyloplace
```

# Running
Run the ```run_phyloplace.sh``` script with the following minimal arguments:
- query_fn - input file
- out_dirn - output directory which will contain: summary file (TSV), FASTA of processed sequences, wise2 prediction results, SVG visualisations of sequences placed into clades of interest (default: macrolactone_ancestral_clade, hr_pks, everything=ALL)

Optional positional arguments (by default, macrolactone reference files are used and included):
- reference protein phylogeny (Newick format)
- reference protein alignment (FASTA)
- reference unaligned protein sequences (FASTA, same sequences)
- reference HMM to be used by WISE2 (HMMer v2 format required, use hmmbuild/hmmconvert on your reference if needed)

Optional keyword arguments:
- annots_fn - annotations TSV file for classification, see data/reference.annots.tsv for macrolactone example)
- clades_fn - clades of interest in 3 column TSV (name of clade, seqid1, seqid2 where seqids are for leaves(=sequences) from reference tree/alignment - see data/clades.tsv for example). Clade of interest is determined by Most Recent Common Ancestor of the two leaves.
- far_dmnd_fn - name of DIAMOND database for further classification of (mainly) outlier sequences, by default SwissProt 2022/05 is included)
- protein - switch indicating query sequences are protein and WISE2 will not be used
- help - display help message
