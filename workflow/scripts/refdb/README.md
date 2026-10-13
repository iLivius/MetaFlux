# Reference-database preparation scripts

This folder holds the **database-specific transform logic** — one script per
reference that MetaFlux cannot use as downloaded — and the check every taxonomy
reference passes before classification (see [The reference check](#the-reference-check)).
Keeping them here (rather than inline in a Snakemake rule) means there is a single,
obvious place to look when you ask *"how is each reference turned into a
classifier-ready file?"*

The rules that call these live in [`workflow/rules/amplicon/10_refdb.smk`](../../rules/amplicon/10_refdb.smk),
which also carries the full **reference preparation map** (which DB is a plain
download, which needs unpacking, which needs one of the scripts below).

| Script | Turns | Into | Used by (rule) |
|---|---|---|---|
| `silva_prepare_trainset.py` | SILVA 144 DADA2 trainset (7 ranks, Domain to Genus) | the same 7-rank DADA2 trainset with every kingdom written without SILVA's `--other` suffix, and the other ranks without it only under `silva_other: strip`; headers only | `convert_silva_train` |
| `silva_to_sintax.py` | SILVA 16S DADA2 trainset (the ranks of the chosen release: 7 on 144, 6 on 138.2; stops on more ranks than letters) | VSEARCH SINTAX FASTA, one letter per rank from the 16S marker pack (`d,k,p,c,o,f,g` on 144, `d,p,c,o,f,g` on 138.2) | `convert_silva_sintax` |
| `frogs_rpob_to_dada2.py` | FROGS rpoB release FASTA (an accession, `Root` and an `[id: N]` tag per rank in every header; the domain written `k__Bacteria`) | DADA2 assignTaxonomy trainset whose header is the 7-rank lineage, with the domain written `d__Bacteria`; stops on any other first rank; headers only | `convert_rpob_to_dada2` |
| `pr2_utax_domain_first.py` | PR2 18S UTAX FASTA (domain under `k:`, supergroup under `d:`) | VSEARCH SINTAX FASTA with the domain under `d:` and the supergroup under `k:`, the order VSEARCH reads; headers only | `convert_pr2_sintax` |

## The reference check

`check_reference_ranks.py` (rule `check_taxonomy_reference`) changes no file. It
reads the one file the classifier is about to use, downloaded, converted or pinned in
the run config, and checks every header against the rank model of the marker pack:
the number of ranks, the prefixes or SINTAX letters and their order, and the names
allowed in the ranks the pack gives rules for (`reference_values`). When the file
does not match, the run stops before any classification and the message names the
rank and quotes example headers. When it matches, the script writes
`stats/taxonomy_reference_check.json`, which `assign_taxonomy` waits for.

## Conventions

- **In-workflow, from canonical sources.** We fetch each database from its upstream
  home and transform it here, rather than re-hosting pre-baked copies. The tree is
  reproducible from a clean clone and tracks upstream when a URL is bumped.
- **Sequence lines are never altered** — only FASTA *headers* are reshaped, so the
  sequences an ASV is compared against are exactly the upstream ones.
- **One script = one database's quirk.** A transform that is generic (a plain
  DADA2↔SINTAX rank remap) should become a single rank-model-aware converter driven
  by the marker pack's `tax_levels` / `rank_letters` / `rank_prefixes` /
  `prefix_style`, not a copy per marker. `silva_to_sintax.py` is deliberately
  SILVA/16S-specific today (see its docstring); generalise it when a second marker
  needs SINTAX.
- **Rank models stay in the marker packs**, not here. These scripts fix *file
  format*; the number of ranks and the prefix style are marker facts that
  `workflow/markers/*.yaml` owns.
