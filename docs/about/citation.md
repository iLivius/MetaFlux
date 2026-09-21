# Citation and references

## Citing MetaFlux

Publications that use MetaFlux should cite the Zenodo record:

> Antonielli, L. (2026). *MetaFlux: a unified short-read multi-marker amplicon and
> shotgun taxonomic profiling workflow.* Zenodo. https://doi.org/10.5281/zenodo.22876451

That DOI always resolves to the newest release. To cite the exact version you ran, use
its own DOI instead — v2.4.0 is
[10.5281/zenodo.22876452](https://doi.org/10.5281/zenodo.22876452) — and give the version
number in the reference, because results can depend on it.

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

MetaFlux was developed in the context of
[MICROBE — MICRObiome Biobanking (RI) Enabler](https://cordis.europa.eu/project/id/101094353),
a Horizon Europe research-infrastructure project coordinated by AIT (2023–2027) that
develops the methods, and the routes of access, for preserving microbiomes with their
composition and function intact. This project has received funding from the European
Union's Horizon Europe research and innovation programme under grant agreement
No. 101094353.

A large part of MetaFlux — merging the amplicon and shotgun workflows into one, the
marker packs, two full code audits, this documentation, and the phylogeny module with
its validation — was built in working sessions with Claude Code. Anthropic accepted
MetaFlux into their Open Source Program and provided Claude Max for it, generously and
with no strings attached; that support is what made the pace and the scope of this work
possible, and it is acknowledged here with real gratitude. The deal on this side has been
that every line is still read, understood and tested by a human before it goes in — which
is why the code and these pages read the way they do.

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

### Phylogeny module (optional, 16S only)

Cite the aligner and whichever tree backend you actually ran — the backend is recorded in
`7.phylogeny/phylogeny.params.json`.

44. Katoh, K. & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30(4):772–780. https://doi.org/10.1093/molbev/mst010 (the aligner, always used when the module runs)
45. Wong, T. K. F., Ly-Trong, N., et al. & Minh, B. Q. (2026). IQ-TREE 3: phylogenomic inference software using complex evolutionary models. *Molecular Biology and Evolution* 43(5):msag117. https://doi.org/10.1093/molbev/msag117 (default tree backend)
46. Kalyaanamoorthy, S., et al. (2017). ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14:587–589. https://doi.org/10.1038/nmeth.4285 (cite only if you set `model: MFP`)
47. Price, M. N., Dehal, P. S. & Arkin, A. P. (2010). FastTree 2 — approximately maximum-likelihood trees for large alignments. *PLoS ONE* 5(3):e9490. https://doi.org/10.1371/journal.pone.0009490 (cite if `backend: fasttree`)
48. Kozlov, A. M., et al. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics* 35(21):4453–4455. https://doi.org/10.1093/bioinformatics/btz305 (cite if `backend: raxml-ng`)
49. Sukumaran, J. & Holder, M. T. (2010). DendroPy: a Python library for phylogenetic computing. *Bioinformatics* 26(12):1569–1571. https://doi.org/10.1093/bioinformatics/btq228 (tree validation, midpoint rooting and the long-branch QC report)

Supporting the module's design decisions, rather than tools it runs:

50. Tan, G., et al. (2015). Current methods for automated filtering of multiple sequence alignments frequently worsen single-gene phylogenetic inference. *Systematic Biology* 64(5):778–791. https://doi.org/10.1093/sysbio/syv033 (why no alignment masking)
51. Janssen, S., et al. (2018). Phylogenetic placement of exact amplicon sequences improves associations with clinical information. *mSystems* 3(3):e00021-18. https://doi.org/10.1128/msystems.00021-18 (the long-branch UniFrac artifact the QC report exists to catch)
52. Mai, U., Sayyari, E. & Mirarab, S. (2017). Minimum variance rooting of phylogenetic trees and implications for species tree reconstruction. *PLoS ONE* 12(8):e0182238. https://doi.org/10.1371/journal.pone.0182238 (midpoint-rooting instability; why the tree is exported unrooted only)

### Methods behind the phylogeny benchmark

Referenced from the [Phylogeny validation](validation.md) page. Not tools MetaFlux runs — the measures and
models used to check its defaults.

53. Robinson, D. F. & Foulds, L. R. (1981). Comparison of phylogenetic trees. *Mathematical Biosciences* 53:131–147. https://doi.org/10.1016/0025-5564(81)90043-2 (Robinson–Foulds distance)
54. Sokal, R. R. & Rohlf, F. J. (1962). The comparison of dendrograms by objective methods. *Taxon* 11:33–40. https://doi.org/10.2307/1217208 (cophenetic / patristic distance correlation)
55. Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. *Biological Conservation* 61:1–10. https://doi.org/10.1016/0006-3207(92)91201-3 (Faith's PD — a sum of branch lengths)
56. Strimmer, K. & von Haeseler, A. (1997). Likelihood-mapping: a simple method to visualize phylogenetic content of a sequence alignment. *PNAS* 94:6815–6819. https://doi.org/10.1073/pnas.94.13.6815
57. Hoang, D. T., et al. (2018). UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35:518–522. https://doi.org/10.1093/molbev/msx281
58. Guindon, S., et al. (2010). New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. *Systematic Biology* 59:307–321. https://doi.org/10.1093/sysbio/syq010 (SH-aLRT branch test)
59. Schwarz, G. (1978). Estimating the dimension of a model. *Annals of Statistics* 6:461–464. https://doi.org/10.1214/aos/1176344136 (BIC, the criterion ModelFinder ranks by)
60. Tavaré, S. (1986). Some probabilistic and statistical problems in the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences* 17:57–86. (GTR model)
61. Kimura, M. (1981). Estimation of evolutionary distances between homologous nucleotide sequences. *PNAS* 78:454–458. https://doi.org/10.1073/pnas.78.1.454 (three-substitution-type models; TPM3u)
62. Yang, Z. (1994). Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods. *Journal of Molecular Evolution* 39:306–314. https://doi.org/10.1007/BF00160154 (discrete gamma rate heterogeneity, +G4)
63. Yang, Z. (1995). A space-time process model for the evolution of DNA sequences. *Genetics* 139:993–1005. (FreeRate site-rate heterogeneity, +R)
64. Soubrier, J., et al. (2012). The influence of rate heterogeneity among sites on the time dependence of molecular rates. *Molecular Biology and Evolution* 29:3345–3358. https://doi.org/10.1093/molbev/mss140 (FreeRate heterogeneity, +R)
65. Capella-Gutiérrez, S., Silla-Martínez, J. M. & Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics* 25:1972–1973. https://doi.org/10.1093/bioinformatics/btp348 (evaluated as an optional mask; not shipped)
66. Bianchini, G., Zhu, Q., Cicconardi, F. & Moody, E. R. R. (2026). AliFilter. *Molecular Biology and Evolution* 43(4):msag097. https://doi.org/10.1093/molbev/msag097 (evaluated as an optional mask; not shipped — not installable from bioconda)

### Phylogeny module — design decisions and the full-gene check

Papers the [Phylogeny](../amplicon/phylogeny.md) page cites for *why* the module is
16S-only, and the [Phylogeny validation](validation.md) page for the full-gene check.

67. Medlar, A., Aivelo, T. & Löytynoja, A. (2014). Séance: reference-based phylogenetic analysis for 18S rRNA studies. *BMC Evolutionary Biology* 14:235. https://doi.org/10.1186/s12862-014-0235-7 (short-fragment branch lengths are underestimated; the 18S deferral)
68. Tedersoo, L., et al. (2022). Best practices in metabarcoding of fungi: from experimental design to results. *Molecular Ecology* 31:2769–2795. https://doi.org/10.1111/mec.16460 (ITS is not alignable beyond genus level)
69. Poirier, S., et al. (2018). Deciphering intra-species bacterial diversity of meat and seafood spoilage microbiota using gyrB amplicon sequencing: a comparative analysis with 16S rDNA V3-V4 amplicon sequencing. *PLoS ONE* 13:e0204629. https://doi.org/10.1371/journal.pone.0204629 (parE paralogue co-amplification; the gyrB exclusion)
70. Case, R. J., et al. (2007). Use of 16S rRNA and rpoB genes as molecular markers for microbial ecology studies. *Applied and Environmental Microbiology* 73:278–288. https://doi.org/10.1128/AEM.01177-06 (rpoB saturation at all three codon positions; the rpoB exclusion)
71. Apprill, A., McNally, S., Parsons, R. & Weber, L. (2015). Minor revision to V4 region SSU rRNA 806R gene primer greatly increases detection of SAR11 bacterioplankton. *Aquatic Microbial Ecology* 75:129–137. https://doi.org/10.3354/ame01753 (806R-B; used with 515F-Y, ref. 27, for the V4 fragments of the full-gene check)

### Downstream analysis in R (the worked example)

Packages and methods the [Phylogeny](../amplicon/phylogeny.md) page's R example uses or
cites. None is a MetaFlux dependency — they run in your own session.

72. Paradis, E. & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics* 35:526–528. https://doi.org/10.1093/bioinformatics/bty633
73. Schliep, K. P. (2011). phangorn: phylogenetic analysis in R. *Bioinformatics* 27:592–593. https://doi.org/10.1093/bioinformatics/btq706 (midpoint rooting)
74. Smith, D. P. & CMMR. rbiom: read/write, transform and visualize microbiome datasets. R package. https://cmmr.github.io/rbiom/ (UniFrac; compiled implementation)
75. Oksanen, J., et al. vegan: Community Ecology Package. R package. https://cran.r-project.org/package=vegan (PCoA, PERMANOVA, betadisper)
76. Davis, N. M., Proctor, D. M., Holmes, S. P., Relman, D. A. & Callahan, B. J. (2018). Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. *Microbiome* 6:226. https://doi.org/10.1186/s40168-018-0605-2 (decontam)
77. Saary, P., Forslund, K., Bork, P. & Hildebrand, F. (2017). RTK: efficient rarefaction analysis of large datasets. *Bioinformatics* 33:2594–2595. https://doi.org/10.1093/bioinformatics/btx206 (multiple rarefaction of non-phylogenetic indices; takes no tree)
78. Lozupone, C. & Knight, R. (2005). UniFrac: a new phylogenetic method for comparing microbial communities. *Applied and Environmental Microbiology* 71:8228–8235. https://doi.org/10.1128/AEM.71.12.8228-8235.2005 · Lozupone, C., Hamady, M., Kelley, S. T. & Knight, R. (2007). Quantitative and qualitative beta diversity measures lead to different insights. *AEM* 73:1576–1585. https://doi.org/10.1128/AEM.01996-06 (weighted UniFrac)
79. Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique and alternative parameters. *Ecology* 52:577–586. https://doi.org/10.2307/1934145 (the rarefaction expectation)
80. Nipperess, D. A. & Matsen, F. A. (2013). The mean and variance of phylogenetic diversity under rarefaction. *Methods in Ecology and Evolution* 4:603–615. https://doi.org/10.1111/2041-210X.12055 (the closed form for expected Faith's PD used in the example)
81. McMurdie, P. J. & Holmes, S. (2014). Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Computational Biology* 10:e1003531. https://doi.org/10.1371/journal.pcbi.1003531 (the case against rarefying)
82. Schloss, P. D. (2024). Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere* 9:e00354-23. https://doi.org/10.1128/msphere.00354-23 · Schloss, P. D. (2024). Waste not, want not: revisiting the analysis that called into question the practice of rarefaction. *mSphere* 9:e00355-23. https://doi.org/10.1128/msphere.00355-23 (the case for rarefaction proper)
83. Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology* 26:32–46. https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x (PERMANOVA) · Anderson, M. J. (2006). Distance-based tests for homogeneity of multivariate dispersions. *Biometrics* 62:245–253. https://doi.org/10.1111/j.1541-0420.2005.00440.x (betadisper)

84. Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research* 27:209–220. (Mantel test for agreement between distance matrices)
85. Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika* 40:33–51. https://doi.org/10.1007/BF02291478 · Peres-Neto, P. R. & Jackson, D. A. (2001). How well do multivariate data sets match? The advantages of a Procrustean superimposition approach over the Mantel test. *Oecologia* 129:169–178. https://doi.org/10.1007/s004420100720 (Procrustes comparison of ordinations; `vegan::protest`)

## License

MetaFlux is released under the
[MIT License](https://github.com/iLivius/MetaFlux/blob/main/LICENSE). Third-party
tools invoked by the workflow are distributed under their own licenses (a mix of
MIT, BSD, and GPL/LGPL); see each tool's repository.
