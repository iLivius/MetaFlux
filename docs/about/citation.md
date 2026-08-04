# Citation and references

## Citing MetaFlux

Publications that use MetaFlux should cite this repository. A formal release with
a Zenodo DOI is forthcoming:

> Antonielli, L. (2026). *MetaFlux: a unified short-read multi-marker amplicon and
> shotgun taxonomic profiling workflow.* Zenodo. DOI: pending release.

MetaFlux is a wrapper around published tools and reference databases, and those
do the actual work. Cite them too — the list below covers the tools and
databases the workflow can invoke, general-purpose helpers such as pigz and seqtk aside,
so pick the entries matching the mode, marker and databases that were actually
used. Two of them cannot be pinned in advance: the Kraken2/Bracken index
(ref. 23) and any host genome supplied for decontamination (ref. 25) are chosen at
runtime, so cite the exact dated build used.

## Acknowledgements

Developed at the [AIT Austrian Institute of Technology](https://www.ait.ac.at/).
MetaFlux consolidates and modernises methods refined across many amplicon and
metagenomics collaborations, and is part of the **BioFlux** family of workflows.

Portions of this codebase were developed with the assistance of Claude Code.

## References

1. Köster, J. & Rahmann, S. (2012). Snakemake — a scalable bioinformatics workflow engine. *Bioinformatics*.
2. Callahan, B. J., et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*.
3. Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*.
4. Langmead, B. & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*.
5. de Sena Brandine, G. & Smith, A. D. (2019). Falco: high-speed FastQC emulation for fastq files. *F1000Research*.
6. Bengtsson-Palme, J., et al. (2015). Metaxa2: improved identification and taxonomic classification of small and large subunit rRNA in metagenomic data. *Molecular Ecology Resources*.
7. Bengtsson-Palme, J., et al. (2013). ITSx: improved software detection and extraction of ITS1 and ITS2. *Methods in Ecology and Evolution*.
8. Rognes, T., et al. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ*.
9. Quast, C., et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*.
10. Abarenkov, K., et al. (2024). UNITE general FASTA release for eukaryotes. *UNITE Community*.
11. Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*.
12. Bushnell, B. BBTools (BBDuk, BBMap). *DOE Joint Genome Institute*.
13. Wood, D. E., Lu, J. & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*.
14. Lu, J., et al. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*.
15. Lu, J., et al. (2022). Metagenome analysis using the Kraken software suite (KrakenTools). *Nature Protocols*.
16. Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*.
17. Wang, Q., et al. (2007). Naïve Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Applied and Environmental Microbiology*.
18. Edgar, R. C. (2016). SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences. *bioRxiv*. https://doi.org/10.1101/074161
19. Guillou, L., et al. (2013). The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote small sub-unit rRNA sequences with curated taxonomy. *Nucleic Acids Research*. (database v5.1.1 — https://github.com/pr2database/pr2database)
20. Briand, M., Rué, O. & Barret, M. (2025). gyrB database for taxonomic assignment formatted for DADA2 (train_set_gyrB_v6). *Recherche Data Gouv (INRAE Dataverse)*. https://doi.org/10.57745/DD7RZ8
21. FROGS rpoB reference databank (2024). Bacterial rpoB genes from complete and chromosome NCBI RefSeq genomes, build 20240707. *INRAE Toulouse*. https://web-genobioinfo.toulouse.inrae.fr/frogs_databanks/assignation/rpoB/
22. Parks, D. H., et al. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research*. (Genome Taxonomy Database, release r226)
23. Kraken 2 / Bracken pre-built index collection. *Langmead Lab, via the AWS Open Data Sponsorship Program.* https://benlangmead.github.io/aws-indexes/k2 (user-supplied at runtime — cite the exact dated build used, e.g. PlusPF)
24. Sanger, F., et al. (1977). Nucleotide sequence of bacteriophage φX174 DNA. *Nature*. (spike-in removal reference: NCBI RefSeq GCF_000819615.1 / NC_001422.1)
25. Handley, S. A. (2020). Virus+ Sequence Masked Human Reference Genome (hg19). *Zenodo*. https://doi.org/10.5281/zenodo.4116107 (default host reference; cite the actual assembly if a different host genome is configured)
26. Stoeck, T., et al. (2010). Multiple marker parallel tag environmental DNA sequencing reveals a highly complex eukaryotic community in marine anoxic water. *Molecular Ecology*. (18S V4 primers TAReuk454FWD1 / TAReukREV3)
27. Parada, A. E., Needham, D. M. & Fuhrman, J. A. (2016). Every base matters: assessing small subunit rRNA primers for marine microbiomes with mock communities, time series and global field samples. *Environmental Microbiology*. (18S V4 primers 515Y / 926R — V4–V5 in 16S nomenclature)
28. Parfrey, L. W., et al. (2014). Communities of microbial eukaryotes in the mammalian gut within the context of environmental eukaryotic diversity. *Frontiers in Microbiology*. (18S V4 primers 515F / 1119r)
29. Hadziavdic, K., et al. (2014). Characterization of the 18S rRNA gene for designing universal eukaryote specific primers. *PLoS ONE*. (18S V4 primers 566F / 1200R)
30. Amaral-Zettler, L. A., et al. (2009). A method for studying protistan diversity using massively parallel sequencing of V9 hypervariable regions of small-subunit ribosomal RNA genes. *PLoS ONE*. (18S V9 primers Euk1391F / EukBr)
31. Barret, M., et al. (2015). Emergence shapes the structure of the seed microbiota. *Applied and Environmental Microbiology*. (gyrB primers F64 / R353; also used to build the DD7RZ8 reference amplicons)
32. Ogier, J.-C., et al. (2019). rpoB, a promising marker for analyzing the diversity of bacterial communities by amplicon sequencing. *BMC Microbiology*. (rpoB primers Univ_rpoB_deg)
33. Dabdoub, S. M. (2016). kraken-biom: enabling interoperative format conversion for Kraken results. *GitHub*. https://github.com/smdabdoub/kraken-biom
34. Leinonen, R., Sugawara, H. & Shumway, M. (2011). The Sequence Read Archive. *Nucleic Acids Research*. (SRA Toolkit — prefetch, fasterq-dump — and NCBI E-utilities; https://github.com/ncbi/sra-tools)
35. Wright, R. J., Comeau, A. M. & Langille, M. G. I. (2023). From defaults to databases: parameter and database choice dramatically impact the performance of metagenomic taxonomic classification tools. *Microbial Genomics* 9(3). https://doi.org/10.1099/mgen.0.000949 (source of the Kraken2 `confidence: 0.15` default)
36. Nyström-Persson, J., Bapatdhar, N. & Ghosh, S. (2025). Precise and scalable metagenomic profiling with sample-tailored minimizer libraries. *NAR Genomics and Bioinformatics* 7(2), lqaf076. https://doi.org/10.1093/nargab/lqaf076 (CAMI2 benchmarking at `confidence: 0.15`)
37. Breitwieser, F. P. & Salzberg, S. L. (2018). KrakenUniq: confident and fast metagenomics classification using unique k-mer counts. *Genome Biology* 19:198. https://doi.org/10.1186/s13059-018-1568-0 (origin of the unique-k-mer evidence idea behind `shotgun.kmer_evidence`)
38. Pochon, Z., et al. (2023). aMeta: an accurate and memory-efficient ancient metagenomic profiling workflow. *Genome Biology* 24:242. https://doi.org/10.1186/s13059-023-03083-9 (source of the paired 1,000 unique k-mers + 200 reads convention)
39. Oskolkov, N. (2026). Refining filtering criteria of Kraken family of tools for accurate taxonomic profiling of ancient metagenomic data. *Frontiers in Microbiology* 17:1603339. https://doi.org/10.3389/fmicb.2026.1603339 (source of `min_distinct_minimizers: 333` after unit conversion, and of `min_reads: 0`)
40. Ye, S. H., et al. (2019). Benchmarking Metagenomics Tools for Taxonomic Classification. *Cell* 178(4):779–794. https://doi.org/10.1016/j.cell.2019.07.010 (principle that a read-count threshold should scale with sequencing depth — the basis of `bracken.threshold: auto`)
41. Meyer, F., et al. (2022). Critical Assessment of Metagenome Interpretation: the second round of challenges. *Nature Methods* 19:429–440. https://doi.org/10.1038/s41592-022-01431-4 (CAMI II marine, plant-associated and strain-madness datasets used to fit `threshold_alpha`)
42. Meyer, F., et al. (2019). Assessing taxonomic metagenome profilers with OPAL. *Genome Biology* 20:51. https://doi.org/10.1186/s13059-019-1646-y (independent scorer used to cross-check the benchmark)
43. Zymo Research. ZymoBIOMICS Gut Microbiome Standard (D6331), Instruction Manual v1.2.0. https://files.zymoresearch.com/protocols/_d6331_zymobiomics_gut_microbiome_standard.pdf (mock community with known composition used for validation)

## License

MetaFlux is released under the
[MIT License](https://github.com/iLivius/MetaFlux/blob/main/LICENSE). Third-party
tools invoked by the workflow are distributed under their own licenses (a mix of
MIT, BSD, and GPL/LGPL); see each tool's repository.
