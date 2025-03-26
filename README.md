# Extensive diversity and rapid turnover of phage defense systems in cheese-associated bacterial communities

Repository containing scripts and statistical analysis of the project. 

*Publication available [here](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01328-6)

<pre><code>

Large variation, rapid turnover and substantial unexplained diversity of phage defense mechanisms and CRISPR in cheese-associated species and communities

Vincent Somerville1,2,*,+, Thibault Schowing2,3,+, Hélène Chabas4, Remo S. Schmidt2,  Ueli von Ah2, Rémy Bruggmann3, Philipp Engel1
1 University of Lausanne, Switzerland
2 Agroscope, Liebefeld, Switzerland
3 Interfaculty Bioinformatics Unit, University of Bern, Bern, Switzerland
4 Institute for Integrative Biology, ETH Zürich, Switzerland
+ These authors contributed equally. 

* Corresponding authors:
vincent.somerville@unil.ch
Department of Fundamental Microbiology
Batiment Biophore
University of Lausanne
1015 Lausanne
Switzerland

</code></pre>






## 1_scripts

Statistical analysis of the pipelines outputs. The outputs consists of the CRISPR content for genomes and metagenomes. Genomes are not available (yet) for property reasons.
Note that the main datasets will be modified during the analysis (e.g. adding sequences GC content / pairwise comparison data / etc) 

| File | Description |
| ------------- | ------------- |
| 1_CRISPRscope_genomic_analysis_phase2.rmd | Genome results R analysis. The pre-processing is done upstream on the HPC side.  |
| 2.1_CRISPRscope_metagenomic_data.rmd      | Metagenomic results pre-processing before statistical analysis (merging and sanitizing)  |
| 2.2_CRISPRscope_metagenomic_analysis.rmd  | Metagenomic results R analysis. |
| 3_CRISPRscope_broad.Rmd                   | R analysis of genomic and metagenomic results together.  |
| 4_BLAST_results_analysis.Rmd              | (Deprecated: moved to 3_CRISPRscope_broad.Rmd) Analysis of BLAST output   |
	
### OtherScripts

Old versions, backups or special scripts.
 
 

## HPC_Scripts


| Directory | Description |
| ------------- | ------------- |
| 0_data_scripts             | Scripts used to retreive data from NCBI or prepare data for BLAST  |
| 1_CRISPRCasFinder          | Main Genomic pipeline. All scripts are chained to produce the raw genomic output.  |
| 2_CRASS                    | Main Metagenomic pipeline.   |
| 3_BLAST                    | BLAST scripts. Note that the BLAST was eventually performed on a different cluster  |
| 4_HMM                      | Markov Models ran on genomic data to annotate defense mechanisms.  |
| 5_SpacerMapping            | Map the DR-SP-DR sequences to raw reads in order to quantify spacers vs protospacers  |
| 6_FastANI                  | Compute the Average Nucleotide Identity of the genomes  |
| 7_metaphlan                | Scripts to get the phylogenetic diversity/content in metagenomes  |
| 8_PairwiseSpacerComparison | Optimised scripts to compute the pairwise sequence identity between the spacers. Performance: ~130 mio comparisons of ~35bp sequences in 2 days with 30 cpus.  |
| 9_GCcontent                | Compute genomes GC content.  |
| A_Clustering               | Given a fasta file, clusters sequences using cd-hit-est and returns a parsed csv output.  |









